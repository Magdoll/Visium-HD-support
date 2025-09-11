# Visium-HD-support

Last updated: 09/10/2025

This document describes how to analyze Visium HD-Kinnex data generated on PacBio long-read sequencers.

From Kinnex sequencing runs, you should use `segmented.bam` (output from SMRT Link Read Segmentation utility (see the [Read Segmentation report section in the single-cell Iso-Seq SMRT Link troubleshooting guide](https://www.pacb.com/wp-content/uploads/SMRT-Link-Kinnex-single-cell-troubleshooting-guide-v13.1.pdf)) or from [skera](https://skera.how/)( via command line) to go into the following steps. 

The `segmented.bam` file should have reads that have sequence IDs in the format `<movie>/<zmw>/ccs/<start>_<end>`, such as `m84039_250124_023230_s2/151066185/ccs/12245_15157`. If the read has sequence ID that is like ``m84039_250124_023230_s2/151066185/ccs`, then it is a HiFi (CCS) reads and still needs to run through Read Segmentation/skera.

If you are not using Kinnex, then a regular HiFi reads BAM file is the proper input. 

## Calling barcodes and UMIs using 10x software
Follow the instructions described [here](https://github.com/10XGenomics/visium-hd-long-reads) to call the barcodes and UMIs using 10x custom scripts and SpaceRanger. The output should consist of a *.spatial_barcodes.csv.gz that has the following format:

```
read_name,uncorrected_barcode,corrected_barcode,uncorrected_umi,corrected_umi
m84039_250124_023230_s2/151066185/ccs/12245_15157,GTCTGCATCTGCCCTGCATTAATGCATCAG,s_002um_02077_014491,CTGGGACGA,CTGGGACGA
m84039_250124_023230_s2/187635340/ccs/2043_2623,GCAGCTATGCAGGTAGTATCCACGGCATCG,s_002um_00964_02405-1,CAATGCATA,CAATGCATA
m84039_250124_023230_s2/146084225/ccs/2498_3082,GCAGCTATGCAGGTAGTATCCACGGCATCG,s_002um_00964_02405-1,CAATGCATA,CAATGCATA
```

## Run a modified single-cell Iso-Seq pipeline
The following workflow is similar to a typical single-cell Iso-Seq workflow (see [example](https://isoseq.how/umi/examples.html)), but we skip the barcode/UMI calling step as that is handled by SpaceRanger.  

### 2.1	Generate the FLNC reads

The following steps trim away the cDNA primers (you can download it [here](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/REF-10x_primers/)), cut out the portions of the sequence attributable to the UMI/BC (but we won’t use the isoseq tag output), remove the polyA tails, to get the full-length non-concatemer (FLNC) reads.
```
$ lima --isoseq --per-read segmented.bam 10x_3kit_primers.fasta <sample>.demux.bam
$ isoseq tag --design T-10U-33B <sample>.demux.5p--3p.bam <sample>.flt.bam
$ isoseq refine <sample>.flt.bam 10x_3kit_primers.fasta <sample>.fltnc.bam --require-polya
```

If you are not using Kinnex libraries, then substitute `segmented.bam` in the first command with HiFi reads BAM file which usually has the file format of `<movie_name>.hifi_reads.bam`.

### 2.2 Custom script to add UMI and CB tags
The following step uses a custom script to add the UMI/BC information from Space Ranger output.

If the library is generated using Kinnex kits, you must firs trun Read Segmentation to get the `segmented.bam` file and provide it here.
```
$ python add_UMI_BC_to_fltnc.py -c <sample>.spatial_barcodes.csv.gz -s <sample>.segmented.bam -i <sample>.fltnc.bam -o <sample> -m <previous retagged CSV output>
```

This will generate a `*<sample>.retagged.BAM*` file. If you have multiple BAM files, you must supply the output of the previous BAM run output `<sample>.retagged.csv` file into the next run, for example:

```
$ python add_UMI_BC_to_fltnc.py -c sample1.spatial_barcodes.csv.gz -s sample1.segmented.bam -i sample1.fltnc.bam -o sample1
$ python add_UMI_BC_to_fltnc.py -c sample2.spatial_barcodes.csv.gz -s sample2.segmented.bam -i sample2.fltnc.bam -o sample2 -m sample1.retagged.csv
$ python add_UMI_BC_to_fltnc.py -c sample3.spatial_barcodes.csv.gz -s sample3.segmented.bam -i sample3.fltnc.bam -o sample3 -m sample2.retagged.csv
...etc
```

If you are using a non-Kinnex library, you do not need to supply the `-s segmented.bam` option.


### 2.3 Run the remaining portion of single-cell Iso-Seq pipeline
We can now run the dedup process on the retagged BAM file.
```
$ samtools sort -t CB <sample>.retagged.bam -o <sample>.corrected.sorted.bam 
$ isoseq groupdedup <sample>.corrected.sorted.bam dedup.bam 
```

Once the reads are deduplicated, we can proceed with alignment to get the isoform GFF files, followed by pigeon isoform classification, and finally count matrix. You can download the relevant hg38 genome and Gencode annotation [here](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/) (if you already have access to SMRT Link, they would also already exist in your SMRT Link bundle).

```
$ pbmm2 align --preset ISOSEQ --sort --report-json mapping_stats.report.json \
    human_GRCh38_no_alt_analysis_set.fasta \
    dedup.bam \
    mapped.bam

$ isoseq collapse mapped.bam collapsed.gff
$ pigeon prepare collapsed.gff
$ pigeon classify collapsed.sorted.gff gencode.v39.annotation.sorted.gtf \
      human_GRCh38_no_alt_analysis_set.fasta
$ pigeon filter collapsed_classification.txt --isoforms collapsed.sorted.gff
$ pigeon make-seurat --dedup dedup.fasta \
     --group collapsed.group.txt \
     -d count_matrix \
     collapsed_classification.filtered_lite_classification.txt
```

### 2.4 Convert sequence back to spatial barcodes
In the `count_matrix/` output directory there are two files that need to be remapped back to spatial barcodes: `genes_seurat/barcode.tsv` and `isoforms_seurat/barcode.tsv`.
Use the following script to convert:
```
$ python remap_to_spatial_barcodes.py -r <sample>.retagged.mapping.csv -i count_matrix/genes_seurat/barcode.tsv
```



