import sys
import pandas as pd
import collections

print ("INP1: Matrix of enzyme substrate associations")
print ("INP2: Output of assign_to_orthogroups.py (og.sp.range) for characterized enzymes")
print ("INP3: Output of assign_to_orthogroups.py (og.sp.range) for test enzymes")

file1=open(sys.argv[1], 'r')
line1=file1.readline()
dict1={}
while line1:
    tab1=line1.strip().split('\t')
    if line1.startswith('#'):
        classes=tab1[1:]        
    else:        
        enz=tab1[0]; binaries=tab1[1:]

        #Get a list of all classes mapped to each enzyme
        list1=[]
        for i in range(0,len(binaries)):
            v1=binaries[i]
            if v1=='1':
                clas=classes[i].replace(' ','_')
                list1.append(clas)

        #Connect positive classes to enzymes
        if enz not in dict1:
            dict1[enz]=list1
        else:
            print ("Enz repeat: ", enz, ", Exiting...")
            #sys.exit()
    line1=file1.readline()
file1.close()

file1=open(sys.argv[2], 'r') #Characterized enzymes range
line1=file1.readline()
dict2={}
while line1:
    tab1=line1.strip().split('\t')
    if line1.startswith('#'):
        pass
    else:        
        enz=tab1[0]; og=tab1[2]; prange=tab1[-1]
        #Associate enzyme with function
        if enz not in dict1:
            print ("Enz absent: ", enz)
            func='NA'
        else:
            func=','.join(dict1[enz])
        str1='{}|{}'.format(enz,func)

        #Associate OG with all enzymes
        if og not in dict2:
            dict2[og]=[str1]
        else:
            if str1 not in dict2[og]:
                dict2[og].append(str1)
    line1=file1.readline()
file1.close()

file1=open(sys.argv[3], 'r') #Novel enzymes range
out1=open(sys.argv[3]+".func",'w')
out2=open(sys.argv[3]+".func.good",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
out2.write('#python {}\n'.format(' '.join(sys.argv)))
out1.write('#Gene\tOG\tMatchHook\tPhylogeneticRange\tMostCommonClass\tAllCharacterized\n')
out2.write('#Gene\tOG\tMatchHook\tPhylogeneticRange\tMostCommonClass\tAllCharacterized\n')

line1=file1.readline()
while line1:
    tab1=line1.strip().split('\t')
    if line1.startswith('#'):
        pass
    else:        
        enz=tab1[0]; og=tab1[2]; prange=tab1[-1]; hit=tab1[1]; sim=tab1[2]
        hook='{}|{}'.format(hit,sim)
        
        if og in dict2:
            enzs=dict2[og]
            str1=';'.join(enzs); slist=[]
            for enzyme in enzs:
                subs=enzyme.split('|')[1].split(',')
                for sub in subs:
                    slist.append(sub)
            counts=collections.Counter(slist)
            mostFreq=counts.most_common(1)

            mf=[]
            for item in mostFreq:
                val=item[0]
                mf.append(val)
            mline=','.join(mf)
        else:
            str1='NA'; mline='NA'
            

        #Find the most frequent substrate class/classes
        out1.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(enz,og,hook,prange,mline,str1))
        if mline!='NA':
            out2.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(enz,og,hook,prange,mline,str1))
    line1=file1.readline()
file1.close(); out1.close(); out2.close()
print ("Done!")
        
            
