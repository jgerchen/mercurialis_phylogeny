from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import itertools
import numpy
#import sys
import argparse


parser=argparse.ArgumentParser(description="Python script to reduce the number of (allelic) sequences in an alignment containing polyploids.")

parser.add_argument("-a", "--alignment", action="store", required=True, help="Input alignment file in fasta format")
parser.add_argument("-p", "--ploidy_file", action="store", required=True, help="Tab seperated file showing the lineage name and number of sequences to be kept for each sample.")
parser.add_argument("-o", "--out_alignment", action="store", required=True, help="Output alignment file in fasta format")
parser.add_argument("-i", "--info_file", action="store", required=True, help="Output file showing which sequences were removed.")


args=parser.parse_args()
arg_dict=vars(args)


aln1=AlignIO.read(arg_dict["alignemnt"], "fasta")
suppl_output=open(arg_dict["out_alignment"], "w")

ploidy_file=open(arg_dict["ploidy_file"])

ploidy_dict={}
for line in ploidy_file:
    ploidy_dict.update({line.strip().split("\t")[0]:int(line.strip().split("\t")[1])})

#TODO: read ploidy dict form input file
#ploidy_dict={"2n":1, "Hue":1, "4n":2, "Can":2, "Rev":2, "Tom":2, "Ell":3, "6n":3, "Per":1}




aln_dict={}
aln_id_count_dict={}

def count_seq_diffs(d_seq1, d_seq2):
    n_diffs=0
    for d_char in range(len(d_seq1)):
        if d_seq1[d_char]!=d_seq2[d_char]:
            n_diffs+=1
    return n_diffs

def reduce_seq_list_by_one(r_seq_list):
    
    suppl_output.write("Reducing list of %s sequences by one\n" % len(r_seq_list))
    suppl_output.write("Sequences: "+",".join(r_seq_list)+"\n")
    seq_combinations=[i for i in itertools.combinations(r_seq_list, 2)]
    #seq_return=seq_list
    comb_differences=[]
    for seq_comb in seq_combinations:
        comb_differences.append(count_seq_diffs(str(aln1[aln_id_count_dict[seq_comb[0]]].seq), str(aln1[aln_id_count_dict[seq_comb[1]]].seq)))
    ###continue here: get index of  combination with lowest differences and randomly remove one of the sequences, return new seq list
    #print(seq_combinations)
    #print(comb_differences)
    
    suppl_output.write("Combinations; N variants:\n")
    for sl_i in range(len(seq_combinations)):
        suppl_output.write("%s - %s; %s\n" % (seq_combinations[sl_i][0], seq_combinations[sl_i][1], comb_differences[sl_i] ))
        
    min_index=numpy.argmin(comb_differences)
    #print(min_index)
    #n_seq_list=r_seq_list
    #n_seq_list.pop(min_index)
    #print(n_seq_list)
    n_r_seq_list=r_seq_list
    remove_seq=numpy.random.choice(seq_combinations[min_index])
    suppl_output.write("Removing sequence %s\n\n" % remove_seq)
    n_r_seq_list.remove(remove_seq)
    return n_r_seq_list

def reduce_copies(sequence_list, copies_wanted):
    #choose random copy
    while len(sequence_list)>copies_wanted:
            sequence_list=reduce_seq_list_by_one(sequence_list)
    return sequence_list


            #n_seq_diffs=count_seq_diffs(str(aln1[aln_id_count_dict[seq_combinations[seq_comb][0]]].seq), str(aln1[aln_id_count_dict[seq_combinations[seq_comb][1]]].seq))
            #print("%s, %s" % (seq_comb, n_seq_diffs))
###########continue here--->count the number of pairwise differences and rank sequences...

aln_dict_counter=0
for aln_seq in aln1:
    #print(aln_seq.id)
    aln_ind=aln_seq.id.split("_")[0] 
    if aln_ind not in aln_dict:
        aln_dict.update({aln_ind:[aln_seq.id]})
    else:
        aln_dict[aln_ind].append(aln_seq.id)
    aln_id_count_dict.update({aln_seq.id:aln_dict_counter})
    aln_dict_counter+=1

output_aln_list=[]

for aln_dict_ind in aln_dict:
    n_copies=0
    if aln_dict_ind[0:2] in ploidy_dict:
        n_copies=ploidy_dict[aln_dict_ind[0:2]]
    elif aln_dict_ind[0:3] in ploidy_dict:
       # print("%s should have %s copies" % (aln_dict_ind, ploidy_dict[aln_dict_ind[0:3]]))
        n_copies=ploidy_dict[aln_dict_ind[0:3]]
    else:
        print("Unknown ploidy!")
    if len(aln_dict[aln_dict_ind])==n_copies:
        copie_counter=1
        print("%s has the right number of copies" % aln_dict_ind)
        for aln_seq in aln_dict[aln_dict_ind]:
            aln_rec=aln1[aln_id_count_dict[aln_seq]]
            aln_rec_id=aln_rec.id
            aln_rec.id=aln_rec_id.strip().split("_")[0]+"_"+str(copie_counter)
            aln_rec.name=aln_rec.id
            aln_rec.description=aln_rec.id
            output_aln_list.append(aln_rec)
            copie_counter+=1
    else:
        print("%s should have %s copies, but has %s" % (aln_dict_ind, n_copies, len(aln_dict[aln_dict_ind]) ))
        new_copies=reduce_copies(aln_dict[aln_dict_ind], n_copies)
        print("%s now has %s copies" % (aln_dict_ind, len(new_copies)) )
        n_copie_counter=1
        for new_copy in new_copies:
            nc_rec=aln1[aln_id_count_dict[new_copy]]
            nc_rec_id=nc_rec.id
            nc_rec.id=nc_rec_id.strip().split("_")[0]+"_"+str(n_copie_counter)
            nc_rec.name=nc_rec.id
            nc_rec.description=nc_rec.id
            output_aln_list.append(nc_rec)
            n_copie_counter+=1
output_aln=MultipleSeqAlignment(output_aln_list)
#TODO: use arparse input for output alignment
#AlignIO.write(output_aln, "subset_fastas/"+aln_name+".fasta", "fasta")


