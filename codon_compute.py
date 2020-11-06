#!/usr/bin/env python3

import os, gzip, itertools

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    
Gene_num1 = 0
len_1 = 0
dict_1 = {'A':0, 'C':0, 'T':0, 'G':0}

numGene_2 = 0
len_2 = 0
dict_2 = {'A':0, 'C':0, 'T':0, 'G':0}

# create dictionary for codons. 
for first in {'A','T','C','G'}:
    for second in {'A','T','C','G'}:
        for third in {'A','T','C','G'}:
            new_codon = first+second+third
            codon[new_codon] = 0
codon_1 = codon.copy()
codon_2 = codon.copy()

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        Gene_num1 += 1
        len_1 += len(seqstring)
        for bp in seq[1]:
            dict_1[bp] += 1
        for n in range(0, len(seq[1]), 3):
            codon_s[seq[1][n:n+3]] += 1
           
GC_1 = (dict_1['G'] + dict_1['C'])/sum(dict_1.values())

with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        Gene_num2 += 1
        len_2 += len(seqstring)
        for bp in seq[1]:
            dict_2[bp] += 1
        for n in range(0, len(seq[1]), 3):
            codon_2[seq[1][n:n+3]] += 1

GC_2 = (dict_2['G'] + dict_2['C'])/sum(dict_2.values())

print("the total number of genes in Salmonella enterica is", Gene_num1)
print("the total num of genes in M. tuberculosis is", Gene_num2)
print("the total length of genes in Salmonella enterica is", len_1)
print("the the total length of genes in M. tuberculosis is", len_2)
print("GC content in Salmonella is",GC_1*100)
print("GC content in Mycobacterium is", GC_2*100)
print("Total number of codon in Salmonella",len_1/3")
print("total number of codon in Mycobacterium is",len_2/3)
