import sys, readFile
print ("INP1: RAP-MSU_2021-05-10.txt")
print ("INP2: uq_norm_rpk.txt.bahd.csv")

file1=open(sys.argv[1],'r')
line1=file1.readline()
dict1={}; dict2={}
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        id1=tab1[0]; locs=tab1[1].split(',')
        #dict2[id1]=[]
        if tab1[1]=='None':
            pass
        else:
            for loc in locs:
                if loc not in dict1:
                    dict1[loc]=id1
            if id1 not in dict2:
                dict2[id1]=[loc]
            else:
                if loc not in dict2[id1]:
                    dict2[id1].append(loc)
    line1=file1.readline()
file1.close()


file1=open(sys.argv[2],'r')
out1=open(sys.argv[2]+".locid",'w')
line1=file1.readline()
y=0; doned={}; d2=[]; x=0
while line1:
    if line1.startswith('id'):
        out1.write(line1)
    else:
        tab1=line1.strip().split('\t')
        id1=tab1[0].replace('t','g').split('-')[0]        
        if id1 in dict2: #List of BAHDs
            locids=dict2[id1]; locid=locids[0]
            print (id1, locids)
            if locid not in doned:
                out1.write('{}\t{}\n'.format(locid,'\t'.join(tab1[1:])))
                doned[locid]=1
            else:
                str1=doned[locid]+1
                doned[locid]=str1
                out1.write('{}-{}\t{}\n'.format(locid,str1,'\t'.join(tab1[1:])))
            #out2.write('{},{}\n'.format(id1,','.join(tab1[1:])))
            #donelist.append(id1)
            y+=1
        else:
            out1.write('{}\t{}\n'.format(id1,'\t'.join(tab1[1:])))
            print ("Absent: ", id1)
            x+=1
            
    line1=file1.readline()
file1.close(); out1.close()
print ("# of BAHD genes in matrix file with loc ID: ", y)
print ("# of BAHD genes in matrix file without loc ID: ", x)

#for gene in nlist:
#    if gene not in donelist:
        
print ("Done!")

            
