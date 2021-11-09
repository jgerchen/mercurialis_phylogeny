import pysam
import primer3
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import argparse

#parse commandline arguments
parser = argparse.ArgumentParser(description="Python script to find primer pairs in conserved regions across bam files.")
parser.add_argument("-s","--seq_list", action='store', required=True, help="List of contigs in the analysis. Tab separated list of contig, linkage group and position in cM.")
parser.add_argument("-g","--genome", action='store', required=True, help="Genome assembly in fasta format.")
parser.add_argument("-b","--bam_list", action='store', required=True, help="List of bam files in the analysis. Tab separated list of bam file and minimum sequencing depth for the file.")
parser.add_argument("-p","--p3_args", action='store', required=True, help="List of primer3 options. Tab separated list of primer3 option and value.")
parser.add_argument("-o","--out_tsv", action='store', required=True, help="Output file in tsv format.")
parser.add_argument("-u","--out_bed", action='store', required=True, help="Output file in bed format.")


args=parser.parse_args()
arg_dict=vars(args)


def read_seq_list(seq_file, genome_assembly):
    genome_dict=SeqIO.to_dict(SeqIO.parse(genome_assembly ,"fasta"))
    seqs={}
    seq_list=open(seq_file)
    for seq in seq_list:
            s_id=seq.strip().split("\t")[0]
            s_lg=seq.strip().split("\t")[1]
            s_pos=seq.strip().split("\t")[2]
            s_seq=str(genome_dict[s_id].seq)
            seqs.update({s_id:(s_lg, s_pos, s_seq)})
    seq_list.close()
    return seqs
    
def read_p3_args(p3_arg_file):
    p3_args=open(p3_arg_file)
    p3_arg_dict={}
    for arg in p3_args:
            option=arg.strip().split("\t")[0]
            value=arg.strip().split("\t")[1]
            if option=='PRIMER_PRODUCT_SIZE_RANGE':
                    range_list=[[int(v2.split("-")[0]), int(v2.split("-")[1])] for v2 in value.split(";")]
                    p3_arg_dict.update({option:range_list})
            elif "." in option:
                    p3_arg_dict.update({option:float(value)})
            else:
                    p3_arg_dict.update({option:int(value)})
    return p3_arg_dict

def read_bam_list(bam_input_file):
    #{bam_name:(bam_file,min_coverage)}
    bam_files={}
    bam_input=open(bam_input_file)
    for bam_file in bam_input:
            b_filename=bam_file.strip().split("\t")[0]
            assert b_filename[-4:]==".bam","bam input filenames must end with .bam!"
            b_min_cov=int(bam_file.strip().split("\t")[1])
            bam_name=b_filename.split("/")[-1][0:-4]
            bam_files.update({bam_name:(pysam.AlignmentFile(b_filename, "rb"), b_min_cov)})
    return bam_files

def get_pileups(bam_inp, rseq, rseq_id, min_block_cov, min_block_size, p3_g_args, output_file, output_bed, rseq_lg, rseq_cm):
    o_file=open(output_file, "a")
    o_bed=open(output_bed, "a")
    i_seqs=[]
    i_block_seqs=[]
    for bam in bam_inp:
        curr_seq=["N" for i in range(len(rseq))]
        curr_block_seq=[0 for i in range(len(rseq))]
        for bp_col in bam_inp[bam][0].pileup(rseq_id):
            #calculate minimum coverage for accepting a site as polymorphic
            min_pol_cov=(2/bam_inp[bam][1])*bp_col.n
            if bp_col.n>=min_block_cov:            
                curr_block_seq[bp_col.pos]=1
            if bp_col.n>=bam_inp[bam][1]:
                loc_dict={"C":0,"A":0,"G":0,"T":0,"-":0,"s":0}
                for pu_read in bp_col.pileups:
                    if pu_read.is_del:
                        loc_dict["s"]+=1
                    elif pu_read.is_refskip:
                        loc_dict["-"]+=1 
                    else:
                        loc_dict[pu_read.alignment.query_sequence[pu_read.query_position]]+=1
                if loc_dict["-"]<min_pol_cov and loc_dict["s"]<min_pol_cov:
                    curr_nucs=[]
                    for nuc in "CAGT":
                        if loc_dict[nuc]>=min_pol_cov:
                            curr_nucs.append(Seq(nuc))
                    if len(curr_nucs)>0:
                        m_nuc=motifs.create(curr_nucs)
                        curr_seq[bp_col.pos]=str(m_nuc.degenerate_consensus)
        i_seqs.append(curr_seq)
        i_block_seqs.append(curr_block_seq)
    #make blocks
    c_blocks=[]
    in_block=False
    curr_block_begin=0
    curr_block_end=0
    for i_site in range(len(rseq)):
        ones=True
        for i_block_seq in i_block_seqs:
            if i_block_seq[i_site]==0:
                ones=False
        if ones==True and in_block==False:
            in_block=True
            curr_block_begin=i_site
        elif ones==False and in_block==True:
            in_block=False
            curr_block_end=i_site-1
            if curr_block_end-curr_block_begin>=min_block_size:
                c_blocks.append((curr_block_begin, curr_block_end))
    if in_block==True:
        curr_block_end=len(rseq)-1
        if curr_block_end-curr_block_begin>=min_block_size:
            c_blocks.append((curr_block_begin, curr_block_end))
            if curr_block_end-curr_block_begin>=min_block_size:
                c_blocks.append((curr_block_begin, curr_block_end))
    #iterate over blocks, make consensus sequence and call Primer3
    if len(c_blocks)>0:
        IUPAC_set=set(["K", "Y", "W", "S", "R", "M", "B", "D", "H", "V"])
        for c_block in c_blocks:
            ambiguity_sites={}
            p3_seq=""
            c_block_pos=c_block[0]
            local_block_position=0
            while c_block_pos<=c_block[1]:
                c_block_nucs=[]
                for i_seq in i_seqs:
                    c_block_nucs.append(i_seq[c_block_pos])
                block_set=set(c_block_nucs)
                if "N" in block_set:
                    p3_seq=p3_seq+"N"
                elif len(block_set)>1 or block_set.isdisjoint(IUPAC_set)==False:
                    p3_seq=p3_seq+"N"    
                    ambiguity_sites.update({local_block_position:block_set})
                else:
                    p3_seq=p3_seq+c_block_nucs[0] 
                c_block_pos+=1
                local_block_position+=1
            #run primer3
            p3_s_args={'SEQUENCE_ID':rseq_id, 'SEQUENCE_TEMPLATE':p3_seq}
            p3_result_dict=primer3.bindings.designPrimers(p3_s_args, p3_g_args)
            n_primers_found=p3_result_dict['PRIMER_PAIR_NUM_RETURNED']
            for primer_pair in range(n_primers_found):
                print("Primer pair %s out of %s" % (primer_pair , n_primers_found ))
                p_left_pos=p3_result_dict["PRIMER_LEFT_"+str(primer_pair)][0]
                p_left_pos_contig=p_left_pos+c_block[0]
                p_left_seq=p3_result_dict["PRIMER_LEFT_"+str(primer_pair)+"_SEQUENCE"]
                p_left_tm=p3_result_dict["PRIMER_LEFT_"+str(primer_pair)+"_TM"]
                p_right_pos=p3_result_dict["PRIMER_RIGHT_"+str(primer_pair)][0]
                p_right_pos_contig=p_right_pos+c_block[0]
                p_right_seq=p3_result_dict["PRIMER_RIGHT_"+str(primer_pair)+"_SEQUENCE"]
		#get reverse complement of right primer
                p_right_seq_rc=str(Seq(p_right_seq, IUPAC.ambiguous_dna).reverse_complement())
                p_right_tm=p3_result_dict["PRIMER_RIGHT_"+str(primer_pair)+"_TM"]
                p_product_size=p3_result_dict["PRIMER_PAIR_"+str(primer_pair)+"_PRODUCT_SIZE"]
                primer_pair_ok=True
                left_N_pos=p_left_seq.find("N")
                if left_N_pos==-1:
                    print("Found left primer without N")
                    l_final_seq=p_left_seq
                else:
                    left_N_pos_block=p_left_pos+left_N_pos
                    if left_N_pos_block in ambiguity_sites:
                        left_a_nucs=ambiguity_sites[left_N_pos_block]
            #check if there is an ambiguity codon in block set
                        left_amb_set=left_a_nucs.intersection(IUPAC_set)
                        left_base_set=left_a_nucs.difference(IUPAC_set)
                        for l_a_s in list(left_amb_set):
                            l_nucs=resolve_ambiguity(l_a_s)
                            for nuc in l_nucs:
                                left_base_set.add(nuc)
                            #now make consensus ambiguity sequence
                        l_pot_sequences=[]
                        for l_base in left_base_set:
                            l_pot_sequences.append(Seq(p_left_seq[0:left_N_pos]+l_base+p_left_seq[left_N_pos+1:]))
                        l_final_motif=motifs.create(l_pot_sequences)
                        l_final_seq=str(l_final_motif.degenerate_consensus)
                        print("Left consensus sequence: %s" % l_final_seq)
                    else:
                        primer_pair_ok=False
                right_N_pos=p_right_seq_rc.find("N")
                if right_N_pos==-1:
                    print("Found right primer without N")
                    r_final_seq=p_right_seq
                    r_final_seq_rc=p_right_seq_rc
                else:
                    right_N_pos_block=p_right_pos-(len(p_right_seq_rc)-right_N_pos)+1
                    if right_N_pos_block in ambiguity_sites:
                        right_a_nucs=ambiguity_sites[right_N_pos_block]
                        right_amb_set=right_a_nucs.intersection(IUPAC_set)
                        right_base_set=right_a_nucs.difference(IUPAC_set)
                        for r_a_s in list(right_amb_set):
                            r_nucs=resolve_ambiguity(r_a_s)
                            for nuc in r_nucs:
                                right_base_set.add(nuc)
                        r_pot_sequences=[]
                        for r_base in right_base_set:
                            r_pot_sequences.append(Seq(p_right_seq_rc[0:right_N_pos]+r_base+p_right_seq_rc[right_N_pos+1:]))
                        r_final_motif=motifs.create(r_pot_sequences)
                        r_final_seq_rc=r_final_motif.degenerate_consensus
                        r_final_seq=r_final_seq_rc.reverse_complement()
                        print("Right consensus sequence: %s" % r_final_seq)
                    else:
                        primer_pair_ok=False
                if primer_pair_ok==True:
                    print("Accepted")
                    o_bed.write("%s\t%s\t%s\t%s\n" % (rseq_id, p_left_pos_contig, p_left_pos_contig+len(l_final_seq), l_final_seq))
                    o_bed.write("%s\t%s\t%s\t%s\n" % (rseq_id, p_right_pos_contig-len(r_final_seq), p_right_pos_contig, r_final_seq_rc))
                    o_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (rseq_id, rseq_lg, rseq_cm, p_product_size, l_final_seq, p_left_pos_contig, p_left_tm, r_final_seq, p_right_pos_contig, p_right_tm))
                else:
                    print("Not accepted")
    o_bed.close()
    o_file.close()

def resolve_ambiguity(amb_dict_inp):
    amb_dict={"M": set(["A", "C"]),"R": set(["A", "G"]),"W": set(["A", "T"]),"S": set(["C", "G"]),"Y": set(["C", "T"]),"K": set(["G", "T"]),"V": set(["A","C", "G"]), "H": set(["A","C","T"]), "D":set(["A", "G", "T"]), "B":set(["C","G","T"])}
    return(amb_dict[amb_dict_inp])

seqs=read_seq_list(arg_dict["seq_list"], arg_dict["genome"])
bams=read_bam_list(arg_dict["bam_list"])
p3_args=read_p3_args("p3_args.tsv")
o_file_name=arg_dict["out_tsv"]
o_file=open(o_file_name, "w")
o_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("scaffold", "lg", "cm", "amplicon_size", "left_seq", "left_pos", "left_tm", "right_seq", "right_pos", "right_tm"))
o_file.close()
o_bed_name=arg_dict["out_bed"]
for seq in seqs:
    print(seq)
    #if len(c_stretch)>0:
    #    print(c_stretch)
    get_pileups(bams, seqs[seq][2], seq, 3, 1000, p3_args,o_file_name , o_bed_name, seqs[seq][0], seqs[seq][1] )


