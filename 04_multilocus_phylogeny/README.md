# Generate multilocus phylogeny

Python script which generates multilocus phylogenies for alignments containing allopolyploid samples with multiple diverged homeologous sequences. It does this by generating a single tree from an alignment and then concatenating a second alignment to the first one, while choosing the combinations of homeologous sequences based on the tree of the first alignment which has the highest log-likelihood for the concatenated alignment.
This file outputs alignemnts (in phylip format) and trees for each iterative step in individual folders.

This script was written for python3 and requires biopython (https://biopython.org/)

In addition it calles RAxML (https://cme.h-its.org/exelixis/web/software/raxml/) for inference of phylogenetic trees and for calculation of log-likelihoods and the path to the RAxML excutable has to be provided.

It requires the following options:

  -a, --alignment_folder Folder containing all alignments used in the analysis in individual files in fasta format

  -r, --raxml_executable Path to the RAxML excutable

  -s, --species_subset Polyploid samples for which all combinations between homeologous sequences between samples will be tested at each individual run. All remaining samples will be added and phased individually one by one, which is much faster. However, it may be possible that the analysis looses the power to accurately phase individuals if the subset of samples (or the diploid samples) does not include individuals that represent each diverged clade in the dataset.
