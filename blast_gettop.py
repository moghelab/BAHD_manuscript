# show only the first match
import sys
file1=open(sys.argv[1], 'r') #blast file
line1=file1.readline()
m=0; d=0
dict1={}
out1=open(sys.argv[1]+".top", 'w')
#out1.write('#python %s\n'%(' '.join(sys.argv)))
while line1:
    tab1=line1.strip().split('\t')
    gene=tab1[0];match=tab1[1]
    if gene!=match:
        if gene not in dict1:        
            dict1[gene]=1
            out1.write(line1)
            d+=1
        else:
            pass
    m+=1
    if m%1000000==0:
        print ("Lines completed: ", m, "Written to OUT: ", d, "Genes: ", len(dict1.keys()))
    line1=file1.readline()
file1.close()        
print ("Done!")    


