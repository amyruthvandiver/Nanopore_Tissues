#!/usr/bin/env python
# coding: utf-8

# In[29]:


#!conda install nanosim
#check install
get_ipython().system('which read_analysis.py')


# In[1]:


#cd /Users/amyvandiver/Box/Nanopore/Analysis/simulated/new


# In[32]:


#for building profile, give nanosim actual fastqs and genome alignment too-- this is real slow 
get_ipython().system('read_analysis.py genome  -i /Users/amyvandiver/Box/Nanopore/Wanagat_Muscle/M166_08232021allpass.fastq -rg /Users/amyvandiver/Box/Nanopore/Timp_data/chrM_rot.fa -a minimap2 -t 11')
       


# In[3]:


from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


# In[4]:


#create template with 0k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
obj_seq=Seq(chrm_str)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq
outfile=open("chrM_del_0k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with1k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[5000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq
outfile=open("chrM_del_1k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 2k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[6000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_2k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 3k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[7000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_3k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 4k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[8000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_4k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 5k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[9000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_5k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 6k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[10000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_6k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 7k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[11000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_7k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 8k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[12000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_8k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 9k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[13000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_9k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()

#create template with 10k del
hg38_chrm=SeqIO.read("/Users/amyvandiver/Box/Nanopore/TimpLab/Ref/chrM.hg38.fa", "fasta")
chrm_str=str(hg38_chrm.seq)
chrm_str=chrm_str[1547:]+chrm_str[0:1547]
del_seq=chrm_str[0:4000]+chrm_str[14000:]
obj_seq=Seq(del_seq)
hg38_del=hg38_chrm
hg38_del.seq=obj_seq

outfile=open("chrM_del_10k.fa", "w")
SeqIO.write(hg38_del, outfile, "fasta")
outfile.close()


# In[6]:


#use nanosim to simulate fastqs based on profile generated 
get_ipython().system('simulator.py genome -rg chrM_del_1k.fa -n 500 -max 16569 -min 100 -o chrM_del_1k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_2k.fa -n 500 -max 16569 -min 100 -o chrM_del_2k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_3k.fa -n 500 -max 16569 -min 100 -o chrM_del_3k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_4k.fa -n 500 -max 16569 -min 100 -o chrM_del_4k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_5k.fa -n 500 -max 16569 -min 100 -o chrM_del_5k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_6k.fa -n 500 -max 16569 -min 100 -o chrM_del_6k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_7k.fa -n 500 -max 16569 -min 100 -o chrM_del_7k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_8k.fa -n 500 -max 16569 -min 100 -o chrM_del_8k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_9k.fa -n 500 -max 16569 -min 100 -o chrM_del_9k -b guppy --perfect --fastq')
get_ipython().system('simulator.py genome -rg chrM_del_10k.fa -n 500 -max 16569 -min 100 -o chrM_del_10k -b guppy --perfect --fastq')


# In[14]:


#define function to align fastqs
def align(fastq, ref,outname):
    #align
    get_ipython().system('minimap2 -ax map-ont {ref} -t 11 {fastq} --MD| samtools view -uF 4 | samtools sort -o {outname}')
    get_ipython().system('samtools index {outname}')


# In[6]:


#Align simulated fastqs, separate for each size
ref="/Users/amyvandiver/Box/Nanopore/Timp_data/chrM_rot.fa"

fastq="chrM_del_1k_aligned_reads.fastq"
align(fastq,ref,"del_1k.bam")

fastq="chrM_del_2k_aligned_reads.fastq"
align(fastq,ref,"del_2k.bam")

fastq="chrM_del_3k_aligned_reads.fastq"
align(fastq,ref,"del_3k.bam")

fastq="chrM_del_4k_aligned_reads.fastq"
align(fastq,ref,"del_4k.bam")

fastq="chrM_del_5k_aligned_reads.fastq"
align(fastq,ref,"del_5k.bam")

fastq="chrM_del_6k_aligned_reads.fastq"
align(fastq,ref,"del_6k.bam")

fastq="chrM_del_7k_aligned_reads.fastq"
align(fastq,ref,"del_7k.bam")

fastq="chrM_del_8k_aligned_reads.fastq"
align(fastq,ref,"del_8k.bam")

fastq="chrM_del_9k_aligned_reads.fastq"
align(fastq,ref,"del_9k.bam")

fastq="chrM_del_10k_aligned_reads.fastq"
align(fastq,ref,"del_10k.bam")


# In[7]:


#combine all into one file
get_ipython().system('cat *.fastq > "all_sim.fastq"')


# In[8]:


#align combined file
fastq="all_sim.fastq"
align(fastq,ref,"all_sim.bam")


# In[11]:


# make subset for igv pic

get_ipython().system('seqtk sample chrM_del_1k_aligned_reads.fastq 10 > 1k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_2k_aligned_reads.fastq 10 > 2k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_3k_aligned_reads.fastq 10 > 3k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_4k_aligned_reads.fastq 10 > 4k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_5k_aligned_reads.fastq 10 > 5k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_6k_aligned_reads.fastq 10 > 6k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_7k_aligned_reads.fastq 10 > 7k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_8k_aligned_reads.fastq 10 > 8k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_9k_aligned_reads.fastq 10 > 9k_sub.fastq')
get_ipython().system('seqtk sample chrM_del_10k_aligned_reads.fastq 10 > 10k_sub.fastq')


# In[12]:


#combine subsetted
get_ipython().system('cat *sub.fastq > "all_sim_sub.fastq"')


# In[17]:


#align subset for IGV pic
fastq="all_sim_sub.fastq"
align(fastq,ref,"all_sim_sub.bam")


# In[ ]:




