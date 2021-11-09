#Python script to cluster basecalled Oxford Nanopore reads from PCR amplicons and generate consensus sequences.

This script was written for python3 and uses the following python packages:

biopython https://biopython.org/

pysam https://github.com/pysam-developers/pysam

scikit-learn https://scikit-learn.org/stable/

numpy https://numpy.org/

In addition, it calls freebayes as an external job.

https://github.com/freebayes/freebayes

It requires a number of input files, which are given with the following options:

  -b, --bamfile Basecalled Oxford Nanopore reads aligned against the genome assembly

  -g, --genome Genome assembly in fasta format.

In addition, there are several required options to indicate the specific genomic position where reads are found as well as parameters for the analysis 

  -c, --contig Genomic contig to analyze

  -s, --start Start position of aligned reads on contig

  -e, --end End positions of aligned reads on contig

  -v, --var_freq_min Minimum frequence to consider a variant (should be adjusted for higher ploidies)

  -p, --ploidy Sample ploidy

  -r, --reads_max Maximum number of reads considered at a variant. Reads at sites with higher sequencing depth will be downsamples

  -m, --min_reads Minimum number of reads considered at a variant. Sites with less reads will not be considered.

Outputfiles are put in an indicated folder

  -o, --out_folder Output folder                                                                             
