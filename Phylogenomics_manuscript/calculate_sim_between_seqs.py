#!/usr/bin/env python
import sys, os
# from collections import defaultdict
# from Bio.Blast.Applications import NcbiblastpCommandline
# from StringIO import StringIO
# from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import numpy as np

### script to calculate the similarity between random pairs of sequences from a supplied fasta file.

file1=open(sys.argv[1],'r') ## fasta file --> all OG BAHDs
number=int(sys.argv[2])

line1=file1.readline() ## make a dictionary of BAHD sequences 1: OG BAHDs
fasta_dict={}
while line1:
    if line1.startswith('>'):
        tab1=line1.strip().split('>')
        enz=tab1[1].split(' ')
        enz=enz[0]
        line1=file1.readline()
        tab1=line1.strip()
        seq=tab1
        fasta_dict[enz]=seq
    else:
        pass
    line1=file1.readline()
file1.close()
print "Number of sequences after OGs in dictionary: ", len(fasta_dict)
#pident_dict={}
pident_list=[]
out1=open(sys.argv[1]+".pident.tab",'w')
out1.write('seqid1\tseqid2\tpident\n')
c1=0
c2=0
i=1
#for _ in range(number):
while i <= number:
    i+=1
    seqid1=random.choice(list(fasta_dict.keys()))
    seqid2=random.choice(list(fasta_dict.keys()))
    #print seqid1, seqid2
    # Create two sequence files
    seq1 = SeqRecord(Seq(fasta_dict[seqid1]),
                   id=seqid1)
    seq2 = SeqRecord(Seq(fasta_dict[seqid2]),
                   id=seqid2)
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")
    # Run BLAST using OS
    os.system('blastp -query seq1.fasta -subject seq2.fasta -outfmt 6 -out output.blp')
    os.system('python2 ~gdm67/scripts/shinhanscripts/ParseBlast.py -f gettop -blast output.blp')
    file2=open('output.blp.top')
    if os.stat("output.blp.top").st_size == 0:
        i-=1
        c2+=1
    line1=file2.readline()
    while line1:
        tab1=line1.strip().split('\t')
        pident=float(tab1[2])
        pident_list.append(pident)
        out1.write('%s\t%s\t%s\n'%(seqid1,seqid2,pident))
        c1+=1
        line1=file2.readline()
    file2.close()
print "Number of comparison's performed: ", len(pident_list),c1,c2
print "Average identity: ", sum(pident_list)/len(pident_list)
print "50th percentile of identity: ", np.percentile(pident_list,50)
print "75th percentile of identity: ", np.percentile(pident_list,75)
print "90th percentile of identity: ", np.percentile(pident_list,90)
print "95th percentile of identity: ", np.percentile(pident_list,95)
out1.close(),file2.close()
