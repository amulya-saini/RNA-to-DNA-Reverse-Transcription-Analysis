#!/usr/bin/env python
# coding: utf-8

# In[1]:


# ignoring warnings
import warnings
warnings.filterwarnings('ignore')

# importing pandas package as pd
import pandas as pd


# In[2]:


# importing Library
import os

# printing the current directory
datapath = os.getcwd()


# In[3]:


# Joining the required file to the datapath
file1_path = datapath + "\\RNA-sequence1.fna"

# Joining the required file to the datapath
file2_path = datapath + "\\RNA-sequence2.fna"


# In[4]:


# importing SeqIO from Biopython
from Bio import SeqIO

#defining a function to read the fasta files and get a list of sequences
def read_fasta(file):
    
    # creating an empty list to store the sequence
    sequences = []
    
    # opening the fasta file as handle
    with open(file, "rt") as handle:

        # reading the sequence
        for record in SeqIO.parse(handle,"fasta"):
            
            # storing the sequence
            sequences.append(str(record.seq))
            
    # returning the extracted sequence
    return sequences


# In[5]:


# reading the sequence 1 fasta file
seq1 = read_fasta(file1_path)

# converting the sequence from list to string
seq1= ''.join(map(str, seq1))


# In[6]:


# import Seq from Biopython package
from Bio.Seq import Seq

# reverse transcribing the RNA seq
cdna_1 = Seq(seq1).reverse_complement().transcribe()

# converting the U's to T's in cDNA
cdna_1 = str(cdna_1).replace('U', 'T')

# print
print("The cDNA sequence is:", cdna_1)


# In[7]:


# Getting the Forward strand by replacing the U's with Ts in RNA seq
forward_strand1 = str(seq1).replace('U', 'T')

# print
print(" The forward strand after reverse transcription:", forward_strand1)


# In[8]:


# determining the 16 combinations of dinucleotides
di_nuc = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]


# In[9]:


# determining the 64 combinations of trinucleotides
tri_nuc = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                  "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                  "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                  "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]


# In[10]:


# creating an empty dictionary to store the frequencies along the di-nucleotides
freq_2 = {}

# for loop to calculate the frequecies
for i in di_nuc:
    
    # counting the frequencies
    count = seq1.count(i)
    
    # adding them to the keys
    freq_2[i] = count


# In[11]:


# creating an empty dictionary to store the frequencies along the tri-nucleotides
freq_3 = {}

# for loop to calculate the frequecies
for i in tri_nuc:
    
    # counting the frequencies
    count = seq1.count(i)
    
    # adding them to the keys
    freq_3[i] = count


# In[12]:


# creating a dataframe from the dictionary
df1 = pd.DataFrame.from_dict(freq_2, orient='index', columns=['Absolute Frequency'])

# calculating the total number of dinucleotides in the given sequence
total = df1['Absolute Frequency'].sum()

# adding a column to the dataframe and calculating the percentage
df1['Percentage'] = (df1['Absolute Frequency'] / total) * 100

# print
df1


# In[13]:


# creating a dataframe from the dictionary
df2 = pd.DataFrame.from_dict(freq_3, orient='index', columns=['Absolute Frequency'])

# calculating the total number of trinucleotides in the given sequence
total = df2['Absolute Frequency'].sum()

# adding a column to the dataframe and calculating the percentage
df2['Percentage'] = (df2['Absolute Frequency'] / total) * 100

# print
df2


# In[14]:


# reading the sequence 2 fasta file
seq2 = read_fasta(file2_path)

# converting the sequence from list to string
seq2= ''.join(map(str, seq2))

# import Seq from Biopython package
from Bio.Seq import Seq

# reverse transcribing the RNA seq
cdna_2 = Seq(seq2).reverse_complement().transcribe()

# converting the U's to T's in cDNA
cdna_2 = str(cdna_2).replace('U', 'T')

# print
print("The cDNA sequence is:", cdna_2)


# In[15]:


# Getting the Forward strand by replacing the U's with Ts in RNA seq
forward_strand2 = str(seq2).replace('U', 'T')

# print
print(" The forward strand after reverse transcription:", forward_strand2)


# In[16]:


# creating an empty dictionary to store the frequencies along the di-nucleotides
freq_2a = {}

# for loop to calculate the frequecies
for i in di_nuc:
    
    # counting the frequencies
    count = seq2.count(i)
    
    # adding them to the keys
    freq_2a[i] = count


# In[17]:


# creating an empty dictionary to store the frequencies along the tri-nucleotides
freq_3a = {}

# for loop to calculate the frequecies
for i in tri_nuc:
    
    # counting the frequencies
    count = seq2.count(i)
    
    # adding them to the keys
    freq_3a[i] = count


# In[18]:


# creating a dataframe from the dictionary
df3 = pd.DataFrame.from_dict(freq_2a, orient='index', columns=['Absolute Frequency'])

# calculating the total number of dinucleotides in the given sequence
total = df3['Absolute Frequency'].sum()

# adding a column to the dataframe and calculating the percentage
df3['Percentage'] = (df3['Absolute Frequency'] / total) * 100

# print
df3


# In[19]:


# creating a dataframe from the dictionary
df4 = pd.DataFrame.from_dict(freq_3a, orient='index', columns=['Absolute Frequency'])

# calculating the total number of trinucleotides in the given sequence
total = df4['Absolute Frequency'].sum()

# adding a column to the dataframe and calculating the percentage
df4['Percentage'] = (df4['Absolute Frequency'] / total) * 100

# print
df4


# In[31]:


# merging the di-nucleotide dataframes on index 
merged_df = pd.merge(df1, df3, left_index=True, right_index=True)

# calculating the fold difference between the percentage columns
merged_df['Fold difference'] = merged_df['Percentage_x'] / merged_df['Percentage_y']

# dropping rows with NaN values as they are zero frequencies and percentages
merged_df = merged_df.dropna()

# print
merged_df


# In[32]:


# mergining the di-nucleotide dataframes on index 
merged_df1 = pd.merge(df2, df4, left_index=True, right_index=True)

# calculating the fold difference between the percentage columns
merged_df1['Fold difference'] = merged_df1['Percentage_x'] / merged_df1['Percentage_y']

# dropping rows with NaN values as they are zero frequencies and percentages
merged_df1 = merged_df1.dropna()

# print
merged_df1


# In[ ]:




