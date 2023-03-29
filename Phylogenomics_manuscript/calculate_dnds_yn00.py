from __future__ import division
import sys, dnds, readFile, random, os
print ("INP1: Coding sequence file")

dict1={}
dict1=readFile.readFasta(sys.argv[1],dict1)
dlist=list(dict1.keys())
os.system('export PATH=/programs/paml4.9j/bin:$PATH')
mainout=open(sys.argv[1]+".rate2",'w')
mainout.write('#python {}\n'.format(' '.join(sys.argv)))
mainout.write('#Seq1\tSeq2\tSeqLen1\tSeqLen2\tdN\tdS\tdNdS\n')
count=0
for i in range(0,len(dlist)):
    for j in range(i+1,len(dlist)):
        n1=dlist[i]; n2=dlist[j]
        s1=dict1[n1]; s2=dict1[n2]
        r=random.randint(0,10000)

        #Make a new seq file
        f1='tmp{}.fa'.format(r)
        tmp1=open(f1,'w')
        l1=len(s1); l2=len(s2)
        tmp1.write('>{}\n{}\n>{}\n{}\n'.format(n1,s1,n2,s2))
        tmp1.close()

        #Perform PRANK        
        os.system('prank -codon -quiet -d={} -o={}.aln > stdout'. \
                  format(f1,f1))

        #Make yn00 control file
        ct=open('yn00.ctl','w')
        ct.write('\tseqfile = {}.aln.best.fas\n'.format(f1))
        ct.write('\toutfile = {}.aln.best.fas.yn\n'.format(f1))
        ct.write('\tverbose = 0\n')
        ct.write('\ticode = 0\n')
        ct.write('\tweighting = 0\n')
        ct.write('\tcommonf3x4 = 0\n')
        ct.write('*\tndata = 1\n')
        ct.close()

        #Run paml
        os.system('/programs/paml4.9j/bin/yn00')
        dnfile=open('2YN.dN').readlines(); dsfile=open('2YN.dS').readlines()
        dn=dnfile[-2].split()[1]
        ds=dsfile[-2].split()[1]
        if float(ds)!=0:
            dnds=float(dn)/float(ds)
        else:
            dnds=float(0)
        mainout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'. \
                   format(n1,n2,len(s1),len(s2),dn,ds,dnds))
        
        #Clear files
        os.system('rm -f yn tmp* 2Y* rst* rub stdout')
        os.system('rm -rf tmp*')
        
        #sys.exit()
mainout.close()
print ("Rates written to out: ", count)
print ("Done!")
        
        
