# Subset alignments

Python script to reduce the number of (allelic) sequences in an alignment containing polyploids. For each sample that has more sequences than indicated in the ploidy file, it counts pairwise sequence differences between between all pairs of sequences and removes a random sequence from the pair with the lowest number of sequence differences until it has the indicated number of sequences.

This script was written for python3 and uses the following python packages:

biopython https://biopython.org/

numpy https://numpy.org/

It requires a number of input files, which are given with the following options:

  -a, --alignment Multiple sequence alignment in fasta format

  -p, --ploidy_file Tab seperated file showing the lineage name (individual names begin with lineage name followed by an underscore and sample name, e.g. lineage name for sample 6n_h2 is 6n) and number of sequences should remain after subsetting. An example file is located in examples/ploidy.list

Outputfiles are put in an indicated folder

  -o, --out_alignment Subsetted alignment file

  -i, --info_file Output file documenting which sequences were removed. 
