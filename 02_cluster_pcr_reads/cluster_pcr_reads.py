import pysam
import itertools
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from pysam import VariantFile
from sklearn.metrics.pairwise import pairwise_distances
import numpy
#import sys
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import argparse

parser=argparse.ArgumentParser(description="Python script to cluster basecalled Oxford Nanopore reads from PCR amplicons and generate consensus sequences.")


parser.add_argument("-b","--bamfile", action="store", required=True, help="Basecalled Oxford Nanopore reads aligned against the genome assembly")
parser.add_argument("-g","--genome", action="store", required=True, help="Genome assembly in fasta format.")
parser.add_argument("-c","--contig", action="store", required=True, help="Genomic contig to analyze.")
parser.add_argument("-s","--start", action="store", required=True, help="Start position of aligned reads on contig.")
parser.add_argument("-e","--end", action="store", required=True, help="End position of aligned reads on contig.")
parser.add_argument("-v","--var_freq_min", action="store", required=True, help="Minimum frequence to consider a variant (should be adjusted for higher ploidies)")
parser.add_argument("-p","--ploidy", action="store", required=True, help="Sample ploidy")
parser.add_argument("-r","--reads_max", action="store", required=True, help="Maximum number of reads considered at a variant. Reads at sites with higher sequencing depth will be downsampled.")
parser.add_argument("-m","--min_reads", action="store", required=True, help="Minimum number of reads considered at a variant. Sites with less reads will not be considered.")
parser.add_argument("-o","--out_folder", action="store", required=True, help="Output folder")



args=parser.parse_args()
arg_dict=vars(args)

#bam_file="minimapped/13.bam"
bam_file=arg_dict["bamfile"]
#ref_fasta="contig_ref.fa"
ref_fasta=arg_dict["genome"]
output_folder=arg_dict["out_folder"]

#contig="scaffold68206"
contig=arg_dict["contig"]
#amplicon_start=1824
amplicon_start=int(arg_dict["start"])
#amplicon_end=2607
amplicon_end=int(arg_dict["end"])
#min_var_freq=0.15
min_var_freq=float(arg_dict["var_freq_min"])
ploidy=int(arg_dict["ploidy"])
max_reads=int(arg_dict["reads_max"])
min_reads=ploidy*int(arg_dict["min_reads"])

do_sequence_clustering=True

def downsample_reads(ds_max_reads, ds_alignment, ds_contig, ds_start, ds_end):
    ds_reads=[dr.query_name for dr in ds_alignment.fetch(ds_contig, ds_start, ds_end)]
    ds_blacklist={}
    if len(ds_reads)>ds_max_reads:
        print("Number of reads at locus %s larger than %s. Blacklisting %s random reads" % (contig,max_reads,len(ds_reads)-max_reads))
        ds_reads_bl=list(numpy.random.choice(ds_reads, len(ds_reads)-ds_max_reads, replace=False))
        for ds_read_bl in ds_reads_bl:
            ds_blacklist.update({ds_read_bl:0})
    return ds_blacklist

def get_variants(vpileup, read_blacklist):
    read_counter=0
    variants={}
    reads={}
    read_list=[]
    for pileupcolumn in vpileup:
        #check if position is masked
        pos_AGCT={"A":[],"G":[],"C":[],"T":[]}
        pos_del=0
        for pileupread in pileupcolumn.pileups:
            count_read=True
            if pileupread.alignment.query_name not in reads:
                #create hash for reads
                if len(reads)<max_reads and pileupread.alignment.query_name not in read_blacklist:
                    reads.update({pileupread.alignment.query_name:read_counter})
                    read_list.append(pileupread.alignment)
                    read_counter+=1
                else:
                    count_read=False
            if not pileupread.is_del and not pileupread.is_refskip and count_read==True:
                p_nuc=pileupread.alignment.query_sequence[pileupread.query_position]
                pos_AGCT[p_nuc].append(pileupread.alignment.query_name)
            elif pileupread.is_del:
                pos_del+=1
        #determine sequencing depth
        p_cov=0
        for p_nuc_var in pos_AGCT:
            p_cov+=len(pos_AGCT[p_nuc_var])
        #determine if variant is present and there is not a high number of deletions
        if p_cov>0 and pos_del<p_cov*min_var_freq:
            for p_nuc_var in pos_AGCT:
                if len(pos_AGCT[p_nuc_var])/p_cov>=min_var_freq and 1-len(pos_AGCT[p_nuc_var])/p_cov>=min_var_freq:
                    nuc_reads={}
                    for p_read in pos_AGCT[p_nuc_var]:
                        nuc_reads.update({p_read:p_nuc_var})
                    if pileupcolumn.pos not in variants:
                        variants.update({pileupcolumn.pos:nuc_reads})
                    else:
                        variants[pileupcolumn.pos].update(nuc_reads)
    return variants, reads, read_list, read_counter
#only does pileup on site list and reads bam files
def get_vars_site_list(var_site_input,bam_input):
    #read bam file
    vs_file=pysam.AlignmentFile(bam_input)
    vs_read_counter=0
    vs_reads={}
    vs_read_list=[]
    vs_variants={}
    vs_pileup=vs_file.pileup(contig, amplicon_start, amplicon_end)
    for vs_pu_column in vs_pileup:
        if vs_pu_column.pos in var_site_input:
            vs_variants.update({vs_pu_column.pos:{}})
            print("Site: %s, Variants: %s" % (vs_pu_column.pos,var_site_input[vs_pu_column.pos]) )
            for vs_puread in vs_pu_column.pileups:
                if vs_puread.alignment.query_name not in vs_reads:
                    vs_reads.update({vs_puread.alignment.query_name:vs_read_counter})
                    vs_read_list.append(vs_puread.alignment)
                    vs_read_counter+=1
                if not vs_puread.is_del and not vs_puread.is_refskip:
                    vs_read_nuc=vs_puread.alignment.query_sequence[vs_puread.query_position]
                    vs_variants[vs_pu_column.pos].update({vs_puread.alignment.query_name:vs_read_nuc})
                else:
                    vs_variants[vs_pu_column.pos].update({vs_puread.alignment.query_name:"INDEL"})
    print("Get_var_sites; vs_variants: %s, vs_reads: %s, vs_read_list: %s, vs_read_counter: %s " % (len(vs_variants), len(vs_reads), len(vs_read_list), vs_read_counter))
    return vs_variants, vs_reads, vs_read_list, vs_read_counter

def get_distance(vc_matrix ,u, i):
    if vc_matrix[u][i]==[0,0]:
        return 0
    else:
        return vc_matrix[u][i][0]/(vc_matrix[u][i][0]+vc_matrix[u][i][1])

def create_distance_matrix(variants, reads, read_list):
        #create variant count matrix
        var_count_matrix=[[[0,0] for u in range(len(read_list))] for i in range(len(read_list))]
        for site in variants:
            read_combos=itertools.combinations(variants[site] , 2)
            for read_combo in list(read_combos):
                if variants[site][read_combo[0]]==variants[site][read_combo[1]]:
                    #determine their position var_count_matrix
                    var_count_matrix[reads[read_combo[0]]][reads[read_combo[1]]][0]+=1
                    var_count_matrix[reads[read_combo[1]]][reads[read_combo[0]]][0]+=1
                else:
                    var_count_matrix[reads[read_combo[0]]][reads[read_combo[1]]][1]+=1
                    var_count_matrix[reads[read_combo[1]]][reads[read_combo[0]]][1]+=1

        return pairwise_distances(numpy.array([[get_distance(var_count_matrix,u, i) for u in range(len(read_list))] for i in range(len(read_list))]), metric='precomputed')

def kmeans_clustering(dist_matrix, ploidy):
    clusterer = KMeans(n_clusters=ploidy)
    cluster_labels=list(clusterer.fit_predict(dist_matrix))
    print([str(i)+":"+str(cluster_labels.count(i)) for i in set(cluster_labels)])
    return cluster_labels

def write_output(o_read_list, o_cluster_labels, o_clusters, o_files):
    for o_list_read_index in range(len(o_read_list)):
        if o_cluster_labels[o_list_read_index] in o_clusters:
            o_clusters[o_cluster_labels[o_list_read_index]].write(o_read_list[o_list_read_index])
    for close_rf_index in o_clusters:
        o_clusters[close_rf_index].close()
        pysam.sort("-o", o_files[close_rf_index]+"_sorted.bam", o_files[close_rf_index]+".bam" )
        pysam.index(o_files[close_rf_index]+"_sorted.bam")

def make_consensus(fasta_list, f_genome, f_contig, f_begin, f_end):
    f_genome_dict=SeqIO.to_dict(SeqIO.parse(f_genome, format="fasta"))
    f_seqs=[list(f_genome_dict[f_contig]) for i in fasta_list]
    f_dict={}
    for f_index in range(len(fasta_list)):
        pysam.sort("-o", fasta_list[f_index][:-4]+"_sort.bam" ,fasta_list[f_index])

        pysam.index(fasta_list[f_index][:-4]+"_sort.bam")
        f_vcf_file_name=fasta_list[f_index][:-4]+"_final.vcf"
        with open(f_vcf_file_name, "w") as f_vcf:
            freeb_call=subprocess.run(["freebayes","-f",f_genome, "-p","1","--min-alternate-fraction", "0.5", "-r", f_contig+":"+str(f_begin)+"-"+str(f_end), fasta_list[f_index][:-4]+"_sort.bam"], check=True, stdout=f_vcf)
        f_vcf_file=pysam.VariantFile(f_vcf_file_name)
        for rec in f_vcf_file.fetch():
            if len(rec.alts[0])==1:
                f_seqs[f_index][rec.pos-1]=rec.alts[0]
            else:
                if len(rec.ref)==len(rec.alts[0]):
                    for i in range(len(rec.ref)):
                        f_seqs[f_index][rec.pos-1+i]=rec.alts[0][i]
                elif len(rec.ref)>len(rec.alts[0]):
                    for i in range(len(rec.ref)):
                        if i<len(rec.alts[0]):
                            f_seqs[f_index][rec.pos-1+i]=rec.alts[0][i]
                        else:
                            f_seqs[f_index][rec.pos-1+i]="-"
                elif len(rec.ref)<len(rec.alts[0]):
                    for i in range(len(rec.ref)):
                        f_seqs[f_index][rec.pos-1+i]=rec.alts[0][i]
        f_dict.update({fasta_list[f_index][:-4]+"_sort.bam":"".join(f_seqs[f_index])[f_begin-1:f_end-1]})
    return f_dict

a_file=pysam.AlignmentFile(bam_file)
a_pileup=a_file.pileup(contig, amplicon_start, amplicon_end)

blacklist_reads=downsample_reads(max_reads, a_file, contig, amplicon_start, amplicon_end)
initial_variants=get_variants(a_pileup, blacklist_reads)
initial_distance_matrix=create_distance_matrix(initial_variants[0], initial_variants[1], initial_variants[2])

initial_km_cluster=kmeans_clustering(initial_distance_matrix,ploidy)
initial_cl_read_counts=[initial_km_cluster.count(i) for i in range(ploidy)]
initial_clusters={}
initial_cl_files={}
initial_n_reads=initial_variants[3]

vcf_var_sites={}
blacklist_clusters=[]

for cl_r_index in range(ploidy):
    if initial_cl_read_counts[cl_r_index]>initial_n_reads*min_var_freq:
        print("Initial cluster %s accepted: %s reads > %s " % (cl_r_index, initial_cl_read_counts[cl_r_index], initial_n_reads*min_var_freq ))
        ic_file=output_folder+"/"+bam_file.split("/")[-1][:-4]+"_"+contig+"_"+str(cl_r_index)+"_in"
        initial_clusters.update({cl_r_index:pysam.AlignmentFile(ic_file+".bam", "wb", template=a_file)})
        initial_cl_files.update({cl_r_index:ic_file})
    else:
        print("Initial cluster %s rejected: %s reads < %s " % (cl_r_index, initial_cl_read_counts[cl_r_index], initial_n_reads*min_var_freq ))
        blacklist_clusters.append(cl_r_index)
write_output(initial_variants[2], initial_km_cluster, initial_clusters, initial_cl_files)

output_fasta=[]

for cl_r_i in range(ploidy):
    if cl_r_i not in blacklist_clusters:
        if do_sequence_clustering==True:
            print("Cluster %s:" % cl_r_i)
            with open(initial_cl_files[cl_r_i]+".vcf","w") as vcfout:
                freeb_call=subprocess.run(["freebayes","-f",ref_fasta, "-C", str(int(0.1*initial_cl_read_counts[cl_r_i]))  ,"-r", contig+":"+str(amplicon_start)+"-"+str(amplicon_end), initial_cl_files[cl_r_i]+"_sorted.bam"], check=True, stdout=vcfout)
            vcf_inp=VariantFile(initial_cl_files[cl_r_i]+".vcf")
            vcf_var_sites.clear()
            for vcf_var in vcf_inp.fetch():
                vcf_gt=[s['GT'] for s in vcf_var.samples.values()][0]
                if len(set(vcf_gt))>1 and vcf_var.info['TYPE'][0]=="snp":
                    vcf_var_sites.update({vcf_var.pos-1:vcf_var.alleles})
            if len(vcf_var_sites)>0:
                vc_out_a=pysam.AlignmentFile(initial_cl_files[cl_r_i]+"_a.bam","wb",template=a_file)
                vc_out_b=pysam.AlignmentFile(initial_cl_files[cl_r_i]+"_b.bam","wb",template=a_file)
                v_site_list=get_vars_site_list(vcf_var_sites ,initial_cl_files[cl_r_i]+"_sorted.bam")
                if len(v_site_list[0])==1:
                    v_position=list(vcf_var_sites.keys())[0]
                for v_site_read in v_site_list[1]:
                    if v_site_read in v_site_list[1] and v_site_list[0][v_position][v_site_read] in vcf_var_sites[v_position]:
                        if len(vcf_var_sites[v_position])==2:
                            if v_site_list[0][v_position][v_site_read]==vcf_var_sites[v_position][0]:
                                vc_out_a.write(v_site_list[2][v_site_list[1][v_site_read]]) 
                            elif v_site_list[0][v_position][v_site_read]==vcf_var_sites[v_position][1]:
                                vc_out_b.write(v_site_list[2][v_site_list[1][v_site_read]]) 
                        elif len(vcf_var_sites[v_position])==3:
                            if v_site_list[0][v_position][v_site_read]==vcf_var_sites[v_position][1]:
                                vc_out_a.write(v_site_list[2][v_site_list[1][v_site_read]]) 
                            elif v_site_list[0][v_position][v_site_read]==vcf_var_sites[v_position][2]:
                                vc_out_b.write(v_site_list[2][v_site_list[1][v_site_read]]) 
                else:
                    print("2nd Clustering, %s; %s; %s" % (len(v_site_list[0]),len(v_site_list[1]),len(v_site_list[2])))
                    vc_distance_matrix=create_distance_matrix(v_site_list[0], v_site_list[1], v_site_list[2])
                    vc_cluster=kmeans_clustering(vc_distance_matrix, 2)
                for vc_read_clust_index in range(len(vc_cluster)):
                    if vc_cluster[vc_read_clust_index]==0:
                        vc_out_a.write(v_site_list[2][vc_read_clust_index]) 
                    elif vc_cluster[vc_read_clust_index]==1:
                        vc_out_b.write(v_site_list[2][vc_read_clust_index]) 
                output_fasta.append(initial_cl_files[cl_r_i]+"_a.bam")
                output_fasta.append(initial_cl_files[cl_r_i]+"_b.bam")
                vc_out_a.close()
                vc_out_b.close()
        else:
            output_fasta.append(initial_cl_files[cl_r_i]+".bam")
output_dict=make_consensus(output_fasta, ref_fasta, contig, amplicon_start, amplicon_end)
fasta_seqs=[]
file_prefix=bam_file.split("/")[-1][:-4]
seq_bam_out=open(output_folder+"/"+file_prefix+"_"+contig+"seq_bam.lst","w")
output_fasta=open(output_folder+"/"+file_prefix+"_"+contig+".fasta","w")
os_counter=1
for output_i in output_dict:
    if output_dict[output_i] not in fasta_seqs:
        fasta_seqs.append(output_dict[output_i])
        fasta_seq_name=file_prefix+"_"+contig+"_"+str(os_counter)
        seq_bam_out.write(fasta_seq_name+"\t"+output_i+"\n")
        output_fasta.write(">"+fasta_seq_name+"\n"+output_dict[output_i]+"\n")
        os_counter+=1
seq_bam_out.close()
output_fasta.close()
