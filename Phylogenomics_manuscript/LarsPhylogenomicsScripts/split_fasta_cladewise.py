#!/usr/bin/env python
import sys, os
from collections import defaultdict

### Script to split a fasta file into different fasta files based on a clade assignment file.


file1=open(sys.argv[1],'r') ## fasta file
file2=open(sys.argv[2],'r') ## clade assignments
line1=file2.readline()
enz_dict=defaultdict(list)
while line1:
    if line1.startswith('enzyme'):
        pass
    else:
        tab1=line1.strip().split('\t')
        enzyme=tab1[0]
        clade1=tab1[1]
        clade2=list(clade1)
        clade=clade2[0]
        enz_dict[clade].append(enzyme)
    line1=file2.readline()    
print len(enz_dict)
line1=file1.readline()
fasta_dict={}
print enz_dict
while line1:
    if line1.startswith('>'):
        tab1=line1.strip().split('>')
        enz=tab1[1]
        line1=file1.readline()
        tab1=line1.strip()
        seq=tab1
        fasta_dict[enz]=seq
    else:
        pass
    line1=file1.readline()
file1.close; file2.close
## write the created dictionaries to new fasta files, one per clade.
out1=open(sys.argv[1]+'.clade_'+sys.argv[3]+'.fasta', 'w')
out2=open(sys.argv[1]+'.clade_'+sys.argv[3]+'.rest.fasta', 'w')
clade=sys.argv[3]
for i in enz_dict[clade]:
    out1.write('>%s\n%s\n'%(i,fasta_dict[i]))
for j in enz_dict:
    if j != clade:
        for i in enz_dict[j]:
            out2.write('>%s\n%s\n'%(i,fasta_dict[i]))
out1.close()
print "fasta files generated"
print "done"
