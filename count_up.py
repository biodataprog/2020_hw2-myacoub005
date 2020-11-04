#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

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



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")
 
#I need to see the file first to know how to work w/ it. I want to know what the gene, start, end columns look like. After looking they are 2,3,4 respectively. 
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        print(row[2],row[3],row[4])
       
#Cool. Now let's get gene number and gene length.
gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"

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



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

#Now I need to find how many genes and gene lenghts
#because I'll be pulling from the same if statement (if row[2] == "gene") I can run these together.
genecount=0
gene_lengths= []
gff = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
with gzip.open(gff, "rt") as fh:
    gff = csv.reader(fh, delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
    if row[2]=="gene":
        genecount += 1
        lengths = int(row[4])-int(row[3])
        gene_lengths.append(lengths)
coding_length = sum(gene_lengths) #this will add all the gene lengths together, giving us total coding length
print("there are",genecount,"in the E. coli genome")
print("the total coding length of the genome is",coding_length)

#Use the FASTA file to compute the total length of genome (by adding up the length of each sequence in the file)

fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

with gzip.open(fasta,"rt") as f:
    seqs = dict(aspairs(f))
total_genome_len = len(seqs['Chromosome']
                      
percent_coding = .format(100* (lengths/total_genome_len)
                         
