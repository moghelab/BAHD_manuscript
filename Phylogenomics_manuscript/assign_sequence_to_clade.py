#!/usr/bin/env python
import fnmatch
import sys, os


### Script to assign sequences to a clade and generate a file for annotating ITOL trees.
### Needs blast output (seq of interest vs characterized), file that assigns colors to each clade,
### and a file with enzyme - clade assignments. Similarity cutoff and length cutoff can be chosen.


file1=open(sys.argv[1],'r') ## clade - color code file
file2=open(sys.argv[2],'r') ## enzyme - clade file
## folder with files
lengthcut=sys.argv[4] ## min length for the hit
simcut=sys.argv[5] ## min similarity cutoff

line1=file1.readline()
color_dict={}
while line1:
    if line1.startswith('color'):
        pass
    else:
        tab1=line1.strip().split('\t')
        color=tab1[0]
        clade=tab1[1]
        color_dict[clade]=color
    line1=file1.readline()    
file1.close
print len(color_dict)
line1=file2.readline()
enz_dict={}
while line1:
    #print line1
    if line1.startswith('enzyme'):
        pass
    else:
        tab1=line1.strip().split('\t')
        enzyme=tab1[0]
        #print enzyme
        clade=tab1[1]
        #print clade
        enz_dict[enzyme]=clade
    line1=file2.readline()    
#print enz_dict
print len(enz_dict)
## read in multiple files to do the clade assigment and write a new 3 column file with enzyme - clade - color
c1=0
lenght_list1=[]
for file in os.listdir(sys.argv[3]): ## blast top hit file
    print file
    if fnmatch.fnmatch(file, '*.top'):
        out1=open(file+'.colors.new1_test.txt', 'w')
        file3=open(file,'r')
        out1.write('TREE_COLORS\n')
        out1.write('SEPARATOR TAB\n')
        out1.write('DATA\n')
        line1=file3.readline()
        while line1:
            tab1=line1.strip().split('\t')
            hit=tab1[0]
            enz=tab1[1]
            perc=tab1[2]
            lenght=tab1[3]
            clade=enz_dict[enz]
            color=color_dict[clade]
            if float(perc) >= float(simcut):
                #print perc
                if int(lenght) >= int(lengthcut):
                    if enz in enz_dict:
                        out1.write('%s\trange\t%s\t%s\n'%(hit,color,clade))
                        #out1.write('%s\tbranch\t%s\tnormal\t1\n'%(hit,color))
                    else:
                        print "ERROR: Enzyme does not exist"
                else:
                    print "hit too short"
                    print lenght
                    pass
            else:
                pass
            line1=file3.readline()
    else:
        print "ERROR: Wrong file name"
out1.close()
