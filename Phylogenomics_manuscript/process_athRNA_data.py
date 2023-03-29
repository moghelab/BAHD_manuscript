import sys
import pandas as pd
print ("INP1: result_table.csv")
file1=open(sys.argv[1],'r')
out1=open(sys.argv[1]+".expr",'w')
line1=file1.readline()
m=0; donelist=[]
while line1:
    if line1.startswith('bad'):
        pass
    elif line1.startswith('Sample'):
        sp=line1.strip().split(',')
        genes=sp[2:66]
        out1.write('SampleName,{}\n'.format(','.join(genes)))
    else:
        sp=line1.strip().split(',')
        #sample='{}_{}'.format(sp[0],sp[1]).replace(' ','_')
        #sample2=sample.replace(', ','_')
        sample=sp[0]
        expr=sp[2:66]
        flag=0
        try:
            for item in expr:
                v1=float(item)
        except:
            flag=1

        if flag==0:
            if sample not in donelist:
                out1.write('{},{}\n'.format(sample,','.join(expr)))        
                m+=1; donelist.append(sample)
        
    line1=file1.readline()
file1.close(); out1.close()
print ("# of genes written: ", m)

f1=pd.read_csv(sys.argv[1]+".expr",sep="\t",header=1,index_col=1)
f2=f1.transpose()
f2.to_csv(sys.argv[1]+".expr.tr",sep="\t")

'''
file1=open(sys.argv[1]+".expr",'r')
line1=file1.readline()
m=0
while line1:
    tab1=line1.strip().split(',')
    m+=1
    if 2390<=m<=2400:
        print (m, len(tab1), tab1)
    if len(tab1)<65:
        print (len(tab1), line1)
    line1=file1.readline()
file1.close()
'''
print ("Done1")
