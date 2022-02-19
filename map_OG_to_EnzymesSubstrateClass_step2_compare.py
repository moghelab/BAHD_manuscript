import sys

print ("INP1: ITAG func file")
print ("INP2: Output of map_OG_step1 (func)")
print ("INP3: Index in INP1 where function is located")
file1=open(sys.argv[1], 'r')
index=int(sys.argv[3])
line1=file1.readline()
dict1={}
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[0]; func=tab1[index]
        if g1 not in dict1:
            dict1[g1]=func
        else:
            print ("Gene repeat 1: ", g1)
    line1=file1.readline()
file1.close()

file1=open(sys.argv[2], 'r')
out1=open(sys.argv[2]+".annotated",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
line1=file1.readline()
dict2={}
while line1:
    if line1.startswith('#'):
        if line1.startswith('#Gene'):
            out1.write('{}\t{}\t{}\n'.format(line1.strip(), 'Gene', 'OriginalAnnot'))
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[0]#.split('.')[0]
        func=tab1[4]; rangex=tab1[3]
        if g1 not in dict2:
            dict2[g1]=f'{func}|{rangex}'
        else:
            print ("Gene repeat 2: ", g1)

        if g1 in dict1:
            ofunc=dict1[g1]
            out1.write('{}\t{}\t{}\n'.format(line1.strip(),g1,ofunc))
        else:
            out1.write('{}\t{}\t{}\n'.format(line1.strip(),g1,'NoPreviousAnnotation'))
            
    line1=file1.readline()
file1.close(); out1.close()

out1=open(sys.argv[2]+".annotated.abbr",'w')
out2=open(sys.argv[2]+".annotated.abbr.good",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
out2.write('#python {}\n'.format(' '.join(sys.argv)))
for gene in dict1:
    flag=0
    f1=dict1[gene]
    if 'HXXXD' in f1 or 'Acyltransf' in f1 or 'Transferase' in f1:
        flag=1

    if gene in dict2:
        f2=dict2[gene]        
    else:
        f2='No_function'
        
    out1.write('{}\t{}\t{}\n'.format(gene,f1,f2))
    if flag==0:
        out2.write('{}\t{}\t{}\n'.format(gene,f1,f2))
        
out1.close(); out2.close()
print ("Done!")
        
