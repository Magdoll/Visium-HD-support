import os, sys
from csv import DictReader
from argparse import ArgumentParser

def remap_barcodes(retagged_csv, input_csv):
    reader = DictReader(open(retagged_csv),delimiter=',')
    d = {}
    for r in reader:
        d[r['CB_tag']+'-1'] = r['spatial_barcode']
    
    f = open(input_csv+'.tmp', 'w')
    for line in open(input_csv):
        f.write(d[line.strip()] + '\n')
    f.close()

    os.rename(input_csv, input_csv+'.old')
    os.rename(input_csv+'.tmp', input_csv)
    print("New barcode CSV written to {0}. Old barcode CSV renamed as {0}.old".format(input_csv))

if __name__ == "__main__":
    parser = ArgumentParser("Remapping sequence back to spatial barcodes")
    parser.add_argument("-r", "--retagged_csv", help="*.retagged.csv file")
    parser.add_argument("-i", "--input_csv", help="Input barcode file name that you want to overwrite")
    args = parser.parse_args()

    remap_barcodes(args.retagged_csv, args.input_csv)
