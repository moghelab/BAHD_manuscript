#!/usr/bin/env python
import sys, os
from collections import defaultdict

### Script that splits a fasta depeding on the assigned orthologous group 



file1=open(sys.argv[1],'r') ## fasta file --> all OG BAHDs
file2=open(sys.argv[2],'r') ## clade assignments
file3=open(sys.argv[3],'r') ## OG assignments
file4=open(sys.argv[4],'r') ## fasta file --> all biochem BAHDs
file5=open(sys.argv[5],'r') ## all_OGs

line1=file2.readline() ## make dictionary of enzymes in clades
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
file2.close
print "Number of BAHDs with clade dictionary: ", len(enz_dict)

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
file1.close
print "Number of sequences after OGs in dictionary: ", len(fasta_dict)
line1=file4.readline() ## add additional BAHD sequences 2: biochem BAHDs 
while line1:
    if line1.startswith('>'):
        tab1=line1.strip().split('>')
        enz=tab1[1].split(' ')
        enz=enz[0]
        line1=file4.readline()
        tab1=line1.strip()
        seq=tab1
        if enz not in fasta_dict:
            fasta_dict[enz]=seq
        else:
            print "double entry: ", enz
    else:
        pass
    line1=file4.readline()
file4.close
print "Number of sequences after biochem BAHDs in dictionary: ", len(fasta_dict)
for i in fasta_dict:
    if len(fasta_dict[i]) == 0:
        print i
        print fasta_dict[i]
line1=file3.readline() ## make dictionary of enzymes in OGs
og_dict=defaultdict(list)
while line1:
    if line1.startswith('enzyme'):
        pass
    else:
        tab1=line1.strip().split('\t')
        enzyme=tab1[0]
        og=tab1[1]
        if og not in og_dict.keys():
            og_dict[og].append(enzyme)
        else:
            og_dict[og].append(enzyme)
    line1=file3.readline()    
file3.close
print "Number of OGs with sequences in dictionary: ", len(og_dict)
#print og_dict
line1=file5.readline() ## make a dictionary with enzymes in OGs
og1_dict=defaultdict(list)
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split(':')
        enzyme_list=tab1[1].strip().split(' ')
        og=tab1[0]
        if og not in og1_dict.keys():
            og1_dict[og]=enzyme_list
        else:
            print "og double entry"
            print og
            print enzyme_list
            print og1_dict[og]
            og1_dict[og]+enzyme_list
    line1=file5.readline()    
file5.close
print "Number of OGs with enzyme entries in dictionary: ", len(og1_dict)
###
### WRITE ROUTINE THAT LOOKS IN THE CLADE DIC, FIND ENZYMES IN THAT CLADE, THEN LOOKS INTO OG DICT AND FINDS OG WITH THAT ENZYME,AFTERWARDS OUTPUT ALL SEQS OF OGS ASSOCIATED WITH
### THAT CLADE
###
c1=0;c2=0;c3=0;c4=0
for i in og1_dict:
    c1 += len(og1_dict[i])
for i in og_dict:
    c2 += len(og_dict[i])
for i in og_dict:
    if i in og1_dict:
        for j in og_dict[i]:
            if j not in og1_dict[i]:
                og1_dict[i].append(j)  
for i in og1_dict:
    c3 += len(og1_dict[i])
for i in og_dict:
    c4 += len(og_dict[i])
print "number of values in og1_dict before: ", c1
print "number of values in og_dict before: ", c2      
print "number of values in og1_dict afer: ", c3
c5=c3-c1
print "Number of enzymes added to multiple OGs ", c5
clade_oglist=[]
i_list=[];j_list=[]
c6=0
for clade in enz_dict:
    print clade
    out1=open(sys.argv[1]+'.'+clade+'_og.fasta', 'w')
    out1.close()
    out2=open(sys.argv[1]+'.'+clade+'_rest.fasta', 'w')
    out2.close()
    clade_oglist=[]
    c=0
    for enzyme in enz_dict[clade]:
        for og in og1_dict:
            if enzyme in og1_dict[og]:
                if og not in clade_oglist:
                    clade_oglist.append(og)
                else:
                    pass
    print clade_oglist
    for og in clade_oglist:
        out1=open(sys.argv[1]+'.'+clade+'_og.fasta', 'a')
        for i in og1_dict[og]:
            if i not in i_list:
                if i in fasta_dict:
                    i_list.append(i)
                    out1.write('>%s\n%s\n'%(i,fasta_dict[i]))
                    c+=1
                else:
                    idd=i.strip().split('_',1)[1]
                    i_list.append(idd)
                    out1.write('>%s\n%s\n'%(idd,fasta_dict[idd]))
                    c+=1
            else:
                if i in i_list:
                    pass
    for og in og1_dict.keys():
        if og not in clade_oglist:
            out2=open(sys.argv[1]+'.'+clade+'_rest.fasta', 'a')
            for i in og1_dict[og]:
                if i in fasta_dict:
                    j_list.append(i)
                    out2.write('>%s\n%s\n'%(i,fasta_dict[i]))
                    c+=1
                else:
                    idd=i.strip().split('_',1)[1]
                    j_list.append(idd)
                    out2.write('>%s\n%s\n'%(idd,fasta_dict[idd]))
                    c+=1
            else:
                pass
    print "Number of sequences written to fasta: ", c
    out1.close();out2.close()
print "done"
