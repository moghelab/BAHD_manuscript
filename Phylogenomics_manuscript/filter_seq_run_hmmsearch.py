from __future__ import division
import sys, readFile, os
print ("INP1: File containing list of files")
print ("INP2: HMM model")

'''
Read all seqs into dict
Count how many seqs less than 100aa
Count how many seqs do not start with M
Add the two
If >20%, fail the dataset
'''

file1=open(sys.argv[1],'r')
out1=open(sys.argv[1]+".passfail",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
line1=file1.readline()
flist=[]; alist=[]
while line1:
    if line1.startswith('#'):
        out1.write('{}\tNotMeth\tLessEqual100\tTot\tPerc\tPassFail\n'.format(line1.strip()))
    else:
        tab1=line1.strip().split('\t')
        fn=tab1[0]; alist.append(fn)
        dict1={}; fail='Passed'
        dict1=readFile.readFasta(fn,dict1)
        print ("Reading: ", fn)

        #Count seqs in each group
        nseq=len(dict1.keys()); nom=0; small=0; flag=0
        for name in dict1:
            seq=dict1[name]
            if seq.startswith('M')==False:
                nom+=1
            elif len(seq)<=100:
                small+=1
            else:
                pass
        #Get total
        flag=nom+small
        perc=(flag/nseq)*100

        #Is it <20%
        if perc>=20:
            fail='Failed'
            flist.append(fn)
            
        out1.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(line1.strip(),nom,small,nseq,perc,fail))
    line1=file1.readline()
file1.close(); out1.close()
print ("# of failed seq datasets: ", len(flist))

#Now perform hmmsearch
file1=open(sys.argv[1]+".passfail",'r')
out1=open(sys.argv[1]+".passfail.count",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
hmm=sys.argv[2]
line1=file1.readline()
while line1:
    if line1.startswith('#'):
        out1.write('{}\tDomainCounts\n'.format(line1.strip()))
    else:
        tab1=line1.strip().split('\t')
        fn=tab1[0]
        
        #Get sequences of all proteins
        sdict={}
        sdict=readFile.readFasta(fn,sdict)
        print sdict
        
        if fn in alist:
            #Find DAOs and count them
            print ("HMM search: ", fn)
            os.system('hmmsearch --noali --cut_tc --cpu 60 --tblout {}.hmm.tab {} {}'. \
                      format(fn,sys.argv[2],fn))

            file2=open(fn+".hmm.tab",'r')
            line2=file2.readline()
            dict2={}
            while line2:
                if line2.startswith('#')==False:
                    name=line2.strip().split()[0]
                    print name
                    dict2[name]=1
                line2=file2.readline()
            file2.close()

            #Write DAO counts to file
            nseqs=len(dict2.keys())
            out1.write('{}\t{}\n'.format(line1.strip(),nseqs))

            #Get DAO sequences
            if fn not in flist:
                out2=open(fn+".hmm.tab.fa",'w')
                for name in dict2:
                    if name in sdict:
                        print name
                        seq=sdict[name]                        
                        out2.write('>{}\n{}\n'.format(name,seq))
                out2.close()
    line1=file1.readline()    
file1.close(); out1.close()

os.system('grep Passed {}.passfail.count > {}.passfail.count.passed'.format(sys.argv[1],sys.argv[1]))
            
print ("Done!")

        
        
        

        
            
