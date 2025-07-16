import re, sys
import pysam
import gzip
from csv import DictReader
from argparse import ArgumentParser

def read_info(spatial_csv_gz, s_bam_file):
    """
    Read the spatial barcodes CSV from SpaceRanger output with format:
    
    read_name,uncorrected_barcode,corrected_barcode,uncorrected_umi,corrected_umi
    m84039_250124_023230_s2/151066185/ccs/12245_15157,GTCTGCATCTGCCCTGCATTAATGCATCAG,s_002um_02077_01449-1,CTGGGACGA,CTGGGACGA
    m84039_250124_023230_s2/187635340/ccs/2043_2623,GCAGCTATGCAGGTAGTATCCACGGCATCG,s_002um_00964_02405-1,CAATGCATA,CAATGCATA

    Also read the segmented.bam to get the `di` tag to trace S-reads

    :returns: spatial_dict, spatial_barcode_map, sread_dict
    """
    spatial_dict = {} # S read name --> CSV dict
    spatial_barcode_map = {} # corrected barcode (not ATCG) --> ATCG-based sequence to use for CB, XC

    # read ahead once to make sure things look right
    rex = re.compile('\S+\/\d+\/ccs\/\d+_\d+')
    with gzip.open(spatial_csv_gz, mode='rt') as h:
        reader = DictReader(h,delimiter=',')
        r = next(reader)
        if not rex.fullmatch(r['read_name']):
            print("Read name in {0} does not have the expected S-read format. Abort!".format(spatial_csv_gz))
            sys.exit(-1)
        if 'corrected_barcode' not in reader.fieldnames:
            print("`corrected_barcode` needs to be a field in {0}. Abort!".format(spatial_csv_gz))
            sys.exit(-1)
        if 'uncorrected_barcode' not in reader.fieldnames:
            print("`uncorrected_barcode` needs to be a field in {0}. Abort!".format(spatial_csv_gz))
            sys.exit(-1)
        if 'corrected_umi' not in reader.fieldnames:
            print("`corrected_umi` needs to be a field in {0}. Abort!".format(spatial_csv_gz))
            sys.exit(-1)

    reader = DictReader(gzip.open(spatial_csv_gz, mode='rt'),delimiter=',')
    for r in reader:
        bc = r['corrected_barcode']
        if bc=='':
            # does not have a corrected barcode, skip!
            #print("Skipping {0} because does not have a corrected barcode.".format(r['read_name']))
            continue
        spatial_dict[r['read_name']] = {'corrected_umi': r['corrected_umi'],
                                        'corrected_barcode': r['corrected_barcode'],
                                        'uncorrected_barcode': r['uncorrected_barcode']}
        if bc not in spatial_barcode_map:
            spatial_barcode_map[bc] = r['uncorrected_barcode'] # we basically use uncorrected bc as a hack for BC/XC
    print("finished reading {0} reads barcode info".format(len(spatial_dict)))

    # read the tagged BAM file to track S-reads by (ZMW, di)
    # where di is the tag that tells you the S-read sub-segment it came from
    # we need this cuz S read names change over time with lima/isoseq refine trimming sequences
    sread_dict = {} # (ZMW,di) --> link to spatial dict
    reader = pysam.AlignmentFile(open(s_bam_file),'rb',check_sq=False)
    for r in reader:
        x = r.qname.split('/')
        zmw = x[0] + '/' + x[1]
        if r.qname not in spatial_dict:
            #print("{0} not in tagged csv, ignore.".format(r.qname))
            continue
        sread_dict[(zmw,dict(r.tags)['di'])] = spatial_dict[r.qname]
    print("finished reading s read info")
    return spatial_dict, spatial_barcode_map, sread_dict


def write_tagged_BAM(in_bam, out_prefix, spatial_dict, spatial_barcode_map, sread_dict):
    """
    Write out a new tagged fltnc with the XM/BC added in
    Input should be fltnc.bam with strand-oriented, cDNA/polyA/UMI-BC trimmed
    The output is equivalent to a barcode corrected BAM that is ready for dedup.
    """
    i = 0
    reader = pysam.AlignmentFile(open(in_bam),'rb',check_sq=False)
    fout = pysam.AlignmentFile(out_prefix+'.retagged.bam', 'wb', header=reader.header)
    for r in reader:
        x = r.qname.split('/')
        tagdict = dict(r.tags)
        _key = x[0]+'/'+x[1], tagdict['di']
        if _key not in sread_dict: 
            print("skipping {0}...".format(_key))
            continue
        i += 1
        newd = r.to_dict()
        new_tags = []
        for tagstr in r.tostring().split('\t')[11:]:
            if tagstr.split(':')[0] not in ['rc','gp','nb','XM','XC','CB','XA']:  new_tags += [tagstr]
        rec = sread_dict[_key]
        new_tags += ['XM:Z:'+rec['corrected_umi'],
                     'CB:Z:'+spatial_barcode_map[rec['corrected_barcode']],
                     'XC:Z:'+rec['uncorrected_barcode'], 
                     'rc:i:1',
                     'gp:i:1',
                     'nb:i:0',
                     'XA:Z:XM-CB']
        newd['tags'] = new_tags
        x = pysam.AlignedSegment.from_dict(newd, r.header)
        fout.write(x)
        #if i > 1000000:
        #    print("HACK BREAK, DELETE THIS LATER")
        #    break
    fout.close()
    
    f = open(out_prefix+'.retagged.mapping.csv', 'w')
    f.write("CB_tag,spatial_barcode\n")
    for spatial_barcode,CB_tag in spatial_barcode_map.items():
        f.write(CB_tag+','+spatial_barcode+'\n')
    f.close()

if __name__ == "__main__":
    parser = ArgumentParser("Retagging BAM file using SpaceRanger barcode information")
    parser.add_argument('-c', '--barcode_csv', help="SpaceRanger barcode CSV file, gzipped")
    parser.add_argument('-s', '--segmented', help="Segmented BAM file")
    parser.add_argument('-i', '--input', help="Input FLTNC BAM file")
    parser.add_argument('-o', '--output_prefix', help="Output prefix")
    
    args = parser.parse_args()
    
    spatial_dict, spatial_barcode_map, sread_dict = read_info(args.barcode_csv, args.segmented)
    write_tagged_BAM(args.input, args.output_prefix, spatial_dict, spatial_barcode_map, sread_dict)
