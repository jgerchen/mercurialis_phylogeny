from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import itertools
import subprocess
import os
from shutil import rmtree
import copy
import re
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
import argparse

parser=argparse.ArgumentParser(description="Python script which generates multilocus phylogenies for alignments containing allopolyploid samples with multiple diverged homeologous sequences.")

parser.add_argument("-a", "--alignment_folder", action="store", required=True, help="Folder which contains all alignments to be used in fasta format.")
parser.add_argument("-r", "--raxml_executable", action="store", required=True, help="Path to the RAxML executable which is called for inferring log-likelihoods and phylogenetic trees.")
#parser.add_argument("-w", "--workdir_path", action="store", required=True, help="Full path to the working directory (required by RAxML).")
parser.add_argument("-s", "--species_subset", action="store", required=True, help="File containing list of samples used for initial analysis with all combinations between samples and homeologs.")

args=parser.parse_args()
arg_dict=vars(args)

#wdir_path=arg_dict["workdir_path"]
wdir_path=os.getcwd()

#check if temp folder exists
if not os.path.exists(wdir_path+"/temp"):
    os.mkdir(wdir_path+"/temp")
else:
    rmtree(wdir_path+"/temp")
    os.mkdir(wdir_path+"/temp")

raxml_path=arg_dict["raxml_executable"]

def calculate_likelihood(cl_alignment, cl_tree):
    raxml_sp=subprocess.run([raxml_path, "-f", "n", "-z", cl_tree, "-s", cl_alignment, "-m", "GTRGAMMA", "-n", "temp", "--silent", "-w", wdir_path+"/temp"], capture_output=True, text=True)
    likelihood=0
    raxml_out=raxml_sp.stdout
    raxml_err=raxml_sp.stderr
    #print(raxml_out)
    likelihood=re.findall("0 (\-[0-9]+\.[0-9]+)",raxml_out)
    #print(likelihood)
    #for raxml_out_line in raxml_out:
        #print(raxml_out_line)
        #if raxml_out_line[0]=="0":
    #likelihood=float(raxml_out_line.strip().split(" ")[1])
    os.remove(wdir_path+"/temp/RAxML_result.temp")
    os.remove(wdir_path+"/temp/RAxML_info.temp")
    os.remove(cl_alignment)
    #os.remove(cl_alignment+".reduced")
    return float(likelihood[0])


def make_permutations(per_samples):
    polyploid_permutations=[]
    perm_sample_list=[]
    #if len(per_samples)>1:
    for per_sample in per_samples:
        if per_samples[per_sample]>1:
            polyploid_permutations.append([u for u in itertools.permutations([per_sample+"_"+str(i) for i in range(1, per_samples[per_sample]+1)], per_samples[per_sample])])
            for i in range(1,per_samples[per_sample]+1):
                perm_sample_list.append(per_sample+"_"+str(i))  
        polyploid_combinations=itertools.product(*polyploid_permutations)
    #else:
    #    per_sample=list(per_samples.keys())[0]
    #    polyploid_combinations=[u for u in itertools.permutations([per_sample+"_"+str(i) for i in range(1, per_samples[per_sample]+1)], per_samples[per_sample])]
    #    for i in range(1,per_samples[per_sample]+1):
    #        perm_sample_list.append(per_sample+"_"+str(i))  
    return perm_sample_list, polyploid_combinations

def subset_alignment(per_samples, s_alignment):
    a_ids={s_alignment[i].id:i for i in range(len(s_alignment))}
    per_sequences={}
    for per_sample in per_samples:
        for ps_i in range(1,per_samples[per_sample]+1):
            per_sequences.update({per_sample+"_"+str(ps_i):0})    
    #make subset of alignments, if there are missing samples (sample from per_samples not in alignment) replace them with Ns
    sub_alignment=MultipleSeqAlignment([])
    for per_sequence in per_sequences:
        if per_sequence in a_ids:
            sub_alignment.append(s_alignment[a_ids[per_sequence]])
        else:
            n_seq=Seq("-"*len(s_alignment[0].seq))
            sub_alignment.append(SeqRecord(n_seq,id=per_sequence))
    return sub_alignment

def add_back_sequences_to_alignment(sub_alignment, full_alignment, add_sequences):
    full_alignment_ids={full_alignment[i].id:i for i in range(len(full_alignment))}
    out_alignment=copy.deepcopy(sub_alignment)
    for add_sequence in add_sequences:
        if add_sequence in full_alignment_ids:
            a_seq=full_alignment[full_alignment_ids[add_sequence]]
        else:
            a_seq=SeqRecord(Seq("N"*len(full_alignment[0])), id=add_sequence)
        out_alignment.append(a_seq)    
    return out_alignment


def subset_tree(st_samples, st_tree, st_tree_file):
    mod_tree=copy.deepcopy(st_tree)
    tree_tips=mod_tree.get_terminals()
    for t_tip in tree_tips:
        if t_tip.name.split("_")[0] not in st_samples:
            mod_tree.prune(t_tip)
    Phylo.write(mod_tree, st_tree_file, "newick")

def cat_permute_alignments(cp_alignment1, cp_alignment2,cp_original, cp_combinations, cp_tree_file):
    assert len(cp_alignment1)==len(cp_alignment2), "Alignments have different numbers of sequences!"
    print("Making permutations...")
    a2_ids={cp_alignment2[i].id:i for i in range(len(cp_alignment2))}
    output_alignments=[]
    comb_counter=0 
    best_likelihood=0
    best_combination=cp_original
    best_alignment=cp_alignment2
    #print(cp_combinations)
    for cp_combination in cp_combinations:
        comb_counter+=1
        if comb_counter%10==0:
            print("Combination %s\n" % comb_counter)
        a2_switch=copy.deepcopy(cp_alignment2)
        cp_combination_chain=[x for x in itertools.chain.from_iterable(cp_combination)]
        for cp_o_i in range(len(cp_original)):
            if cp_original[cp_o_i]!=cp_combination_chain[cp_o_i]:
                a2_switch[a2_ids[cp_original[cp_o_i]]].seq=cp_alignment2[a2_ids[cp_combination_chain[cp_o_i]]].seq
        ###TODO: continue here
        AlignIO.write(cp_alignment1+a2_switch, wdir_path+"/temp/temp.phy", "phylip")  
        curr_likelihood=calculate_likelihood(wdir_path+"/temp/temp.phy", wdir_path+"/temp/temp.nwk")
        if best_likelihood==0 or curr_likelihood>best_likelihood:
            best_likelihood=curr_likelihood
            best_combination=cp_combination
            best_alignment=a2_switch
    print("Analyzed %s combinations. Best likelihood: %s" % (comb_counter, best_likelihood))
    return best_likelihood, best_combination, best_alignment

def phase_two_loci(l1_alignment, l2_alignment, l1_tree, combination_samples, all_samples):
    reduced_permutations=make_permutations(combination_samples)
    print("Made permutations of combined samples...")
    haploid_samples={i:all_samples[i] for i in all_samples if all_samples[i]==1}
    sub_samples={**haploid_samples, **combination_samples}
    sub_alignment_1=subset_alignment(sub_samples, l1_alignment)
    sub_alignment_2=subset_alignment(sub_samples, l2_alignment)
    sub_tree=subset_tree(sub_samples, l1_tree, wdir_path+"/temp/temp.nwk")
    print("Phasing combined samples...")
    per_alignment=cat_permute_alignments(sub_alignment_1, sub_alignment_2,reduced_permutations[0] ,reduced_permutations[1], wdir_path+"/temp/temp.nwk")[2]
    # add missing polyploids and phase them one by one
    remaining_samples=set(all_samples.keys()).difference(set(sub_samples.keys()))
    for remaining_sample in remaining_samples:
        print("Adding sample %s and phasing individually" % remaining_sample)
        sub_samples.update({remaining_sample:all_samples[remaining_sample]})
        rs_permutations=make_permutations({remaining_sample:all_samples[remaining_sample]})
        sub_alignment_1=add_back_sequences_to_alignment(sub_alignment_1,l1_alignment,rs_permutations[0])
        per_alignment=add_back_sequences_to_alignment(per_alignment,l2_alignment,rs_permutations[0])
        # make tree with missing tip added
        rs_tree=subset_tree(sub_samples, l1_tree, wdir_path+"/temp/temp.nwk")
        per_alignment=cat_permute_alignments(sub_alignment_1, per_alignment, rs_permutations[0], rs_permutations[1], wdir_path+"/temp/temp.nwk")[2]
    l1_alignment.sort()
    per_alignment.sort()
    return l1_alignment+per_alignment

def raxml_tree_search(input_alignment, raxml_workdir, raxml_filename):
    #make workdir
    
    if not os.path.exists(raxml_workdir):
        os.mkdir(raxml_workdir)
    else:
        rmtree(raxml_workdir)
        os.mkdir(raxml_workdir)
    AlignIO.write(input_alignment,raxml_workdir+"/"+raxml_filename+".phy","phylip")
    print(wdir_path+raxml_workdir)
    #submit multiprocess raxml job to slurm?
    raxml_sp=subprocess.run([raxml_path, "-f", "a", "-x","123435", "-p", "12345", "-#", "1000", "-s", raxml_workdir+"/"+raxml_filename+".phy", "-m", "GTRGAMMA", "-n", raxml_filename, "--silent", "-w", wdir_path+"/"+raxml_workdir], capture_output=True, text=True)
    tree_result=Phylo.read(raxml_workdir+"/RAxML_bipartitions."+raxml_filename, "newick")
    return tree_result




def analyze_multiple_alignments(alignment_files, comb_samples):
    alignments={}
    samples={}
    combination_samples={line.strip():0 for line in open(comb_samples)}
    for a_file in alignment_files:
        a_name=a_file.split("/")[-1].replace(".fasta", "")
        alignments.update({a_name:AlignIO.read(a_file, "fasta")})
        for a_sample in alignments[a_name]:
            sample_seq=a_sample.id
            sample_name=sample_seq.split("_")[0]
            sample_index=int(sample_seq.split("_")[1])
            if sample_name not in samples:
                samples.update({sample_name:sample_index})
            elif sample_index>samples[sample_name]:
                samples[sample_name]=sample_index
    for combination_sample in combination_samples:
        combination_samples[combination_sample]=samples[combination_sample]
    #sort alignments by number of sequences and then by length
    alignments_sorted=sorted(alignments.keys(), key=lambda x:(len(alignments[x]), len(alignments[x][0].seq)), reverse=True)        
    #make first tree from largest alignment
    tree_counter=0
    phased_aln=alignments[alignments_sorted[0]]
    for sorted_alignment in alignments_sorted:
        tree_directory="Tree_"+str(tree_counter)+"_"+sorted_alignment
        if tree_counter==0:
            #make initial tree
            print("Generating initial tree...")
            prev_tree=raxml_tree_search(alignments[sorted_alignment], tree_directory, str(tree_counter))
            print("Made initial tree...")
        else:
            #otherwise phase sequences first
            print("Phasing alignment combination %s..." % tree_counter)
            phased_aln=phase_two_loci(phased_aln, alignments[sorted_alignment],prev_tree, combination_samples, samples)
            print("Generating Tree %s..." % tree_counter)
            prev_tree=raxml_tree_search(phased_aln,tree_directory, str(tree_counter))
        tree_counter+=1



input_alignments=glob.glob(arg_dict["alignment_folder"]+"/*")
analyze_multiple_alignments(input_alignments ,arg_dict["species_subset"])
