Python script to find primer pairs in conserved regions across bam files.

This script was written for python3 and uses the following python packages:

biopython https://biopython.org/
pysam https://github.com/pysam-developers/pysam
primer3-py https://github.com/libnano/primer3-py

It requires a number of input files, which are given with the following options:

  -s, --seq_list List of contigs to be used in the analysis. Tab separated list of contig, linkage group and position in cM.
	An example file is located in examples/seq.list

  -g, --genome Genome assembly in fasta format.

  -b, --bam_list List of bam files in the analysis. Tab separated list of bam file and minimum sequencing depth for the file.
	An example file is located in examples/bam.list

  -p, --p3_args List of primer3 options. Tab separated list of primer3 option and value.
	An example file is located in examples/p3.list

This script produces two output files, a detailed tsv file and a more reduced bed file, which can be loaded in IGV

  -o, --out_tsv Output file in tsv format.                                                                             
  -u, --out_bed Output file in bed format.  
