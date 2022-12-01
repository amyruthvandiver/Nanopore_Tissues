#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Define reference, set up align function

ref="/Users/amyvandiver/Box/Nanopore/Timp_data/chrM_rot.fa"

def align(fastq, ref,outname):
    #align
    get_ipython().system('minimap2 -ax map-ont {ref} -t 11 {fastq} --MD| samtools view -uF 4 | samtools sort -o {outname}')
    get_ipython().system('samtools index {outname}')


# In[1]:


get_ipython().system('mkdir Wanagat_Muscle')


# In[2]:


cd Wanagat_Muscle/


# In[4]:


import pandas as pd


# In[5]:


#load annotation for muscle samples
anno=pd.read_csv("anno.csv")


# In[19]:


get_ipython().system('mkdir bams')


# In[27]:


# align muscle samples 
for x in anno.Run:
    print(x)
    gzs="/Users/amyvandiver/Box/2021_amy_vandiver/Sequencing_Data/Austin/Muscle_Cas9_Nanopore/cas9_minion_"+x+"/fastq_pass/*.fastq.gz"
    fqs="/Users/amyvandiver/Box/2021_amy_vandiver/Sequencing_Data/Austin/Muscle_Cas9_Nanopore/cas9_minion_"+x+"/fastq_pass/*.fastq"
    file=x+"allpass.fastq"
    get_ipython().system('gunzip {gzs}')
    get_ipython().system('cat {fqs} > {file}')
    get_ipython().system('echo $(( $(wc -l < $file) / 4 )) reads')
    
    bamname="bams/"+x+".rot.bam"
    align(file, ref,bamname)


# In[49]:


# count supplemental and prim alignments in each
for x in anno.Run:
    bamname="bams/"+x+".rot.bam"
    print(x)
    get_ipython().system('samtools view -f 2048 {bamname} -c')
    get_ipython().system('samtools view -F 260 {bamname} -c')


# In[ ]:


#load annotation for brain samples 
anno=pd.read_csv("anno_brain.csv")


# In[ ]:


#align all brain samples 
for x in anno.Run:
    print(x)
    gzs="/Users/amyvandiver/Box/2021_amy_vandiver/Sequencing_Data/Austin/Other_Cas9_Nanopore/cas9_minion_"+x+"/fastq_pass/*.fastq.gz"
    fqs="/Users/amyvandiver/Box/2021_amy_vandiver/Sequencing_Data/Austin/Other_Cas9_Nanopore/cas9_minion_"+x+"/fastq_pass/*.fastq"
    file=x+"allpass.fastq"
    get_ipython().system('gunzip {gzs}')
    get_ipython().system('cat {fqs} > {file}')
    get_ipython().system('echo $(( $(wc -l < $file) / 4 )) reads')
    
    bamname="bams/"+x+".rot.bam"
    align(file, ref,bamname)


# In[ ]:


#count supp and prim alignments for brain samples

for x in anno.Run:
    bamname="bams/"+x+".rot.bam"
    print(x)
    get_ipython().system('samtools view -f 2048 {bamname} -c')
    get_ipython().system('samtools view -F 260 {bamname} -c')

