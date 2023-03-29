import sys, readFile, os
'''
Read in BAHD list
For each BAHD, get coexp genes above z-score of 3 (>99% of the correlations)
Get the pathways of those genes
Do enrichment test COEXP+PWY+   COEXP+PWY-   COEXP-PWY+   COEXP-PWY-
Perform multiple testing correction
Identify enriched pathways for each BAHD gene
'''
print ("INP1: AT BAHD list")
print ("INP2: Entrez-TAIR-PWY correspondence")
print ("INP3: Directory where the ATTED files are located")

#Get a list of all BAHDs
dict1={}
dict1=readFile.readIndex(sys.argv[1], dict1, 0)

#Make a dict of all BAHDs and their ENTREZ values
dict2={}; dict3={}; dict22={}
dict2=readFile.readTwoIndexes(sys.argv[2], dict2, 1, 0)
dict22=readFile.readTwoIndexes(sys.argv[2], dict22, 0, 1)
dict3=readFile.readTwoIndexes(sys.argv[2], dict3, 0, 3)
'''
count=0
for gene in dict3:
    if dict3[gene]!=['NA']:
        count+=1
#gene1=list(dict3.keys())[0]
print (dict3[gene1])
print (count)
sys.exit()
'''
wdir=sys.argv[3]
noentrez=0; dict4={}; dict42={}
##########
def makeStr(gx,list1):
    xstr1='$$'.join(list1)
    xstr2='{}=={}'.format(gx,xstr1)
    return xstr2
##########
print ("Reading ATTED files and getting coexp. and PWY info...")
dlist=list(dict1.keys())#[0:10]
pwydict={}
out1=open(sys.argv[1]+".coexp.pwys",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
for bahd in dlist:
    #Check if BAHD has an Entrez ID
    if bahd in dict2:
        entrezx=dict2[bahd]
        if len(entrezx)==1:
            entrez=entrezx[0]; v1=bahd
            #Get PWYs of the BAHD
            if entrez in dict3:
                pwysb=dict3[entrez]
                str2=makeStr(bahd,pwysb)
                dict4[v1]=[str2]
                for p in pwysb:
                    pwydict[p]=1

            #Get PWYs of coexp genes
            file1=open('{}/{}'.format(wdir, entrez), 'r')
            print ("Reading: ", entrez, bahd)
            line1=file1.readline()
            while line1:
                tab1=line1.strip().split('\t')
                v2=tab1[0]; v3=float(tab1[1])
                if v1!=v2:                
                    #Only add highly co-exp hits with PWY info to final list
                    if v3>=3.0: #(highly co-expressed)
                        #Check if the coexp gene has TAIR ID and PWY info
                        if v2 in dict22 and v2 in dict3:
                            at2=dict22[v2][0]
                            pwys=dict3[v2]                        
                            str2=makeStr(at2,pwys)
                            out1.write('{}\t{}\tCOEXP\t{}\t{}\n'.format(v1,at2,v3,pwys))
                            for p in pwys:
                                pwydict[p]=1
                                
                            #Add gene and PWY info to dict4
                            if v1 not in dict4:
                                dict4[v1]=[str2]
                            else:
                                if str2 not in dict4[v1]:
                                    dict4[v1].append(str2)                        

                    elif -1<=v3<=1.0: #CERTAINLY Not co-expressed (for fisher-exact test)
                        if v2 in dict22 and v2 in dict3:
                            at2=dict22[v2][0]
                            pwys=dict3[v2]
                            str2=makeStr(at2,pwys)                            
                            for p in pwys:
                                pwydict[p]=1
                            
                            #Add gene and PWY info to dict42
                            if v1 not in dict42:
                                dict42[v1]=[str2]
                            else:
                                if str2 not in dict42[v1]:
                                    dict42[v1].append(str2)
                    else:
                        pass
                    
                line1=file1.readline()
            file1.close()
            
    else:
        noentrez+=1
#COEXPRESSED GENES -- dict4[BAHD]=['COEXPGENE==PWY1$$PWY2$$PWY3','COEXPGENE2==PWY3$$PWY4']
#NOT-COEXPRESSED GENES -- dict42[BAHD]=[COEXPGENE==PWY1$$PWY2$$PWY3]
#LIST OF ALL PATHWAYS -- pwydict[PWYNAME]
print ("# of BAHDs with coexp > 3: ", len(dict4.keys()))
#sys.exit()
'''
Do enrichment test COEXP+PWY+   COEXP+PWY-   COEXP-PWY+   COEXP-PWY-
'''
plist=list(pwydict.keys())#.remove("NA")
#print (plist)
glist=dlist
countd={}; totcount=0; doned={}
out2=open(sys.argv[1]+".coexp.fisherinp",'w')
out2.write('#Gene==PWY\tCOEXP+PWY+\tCOEXP+PWY-\tCOEXP-PWY+\tCOEXP-PWY-\n')
out3=open(sys.argv[1]+".coexp.fisherinp.full",'w')
out3.write('#Gene==PWY\tCOEXP+PWY+\tCOEXP+PWY-\tCOEXP-PWY+\tCOEXP-PWY-\n')

print ("Making file for performing enrichment test...")
print ("# of pathways to scan per gene: ", len(plist))
for gene in glist: #all BAHDs
    print ("Analyzing: ", gene)
    for pwy in plist: #list of all pathways
        if pwy!="NA":
            str1='{}=={}'.format(gene,pwy)
            v1d={}; v2d={}; v3d={}; v4d={}            
            
            if gene in dict4: #if it has co-expressed genes
                cgenes=dict4[gene]            
                for g1x in cgenes:
                    g1=g1x.split('==')[0]                
                    if g1!=gene:                    
                        pwys=g1x.split('==')[1].split('$$')                    
                        if pwy in pwys:
                            v1d[g1]=1
                        else:
                            v2d[g1]=1
            
            if gene in dict42: #if it has not-co-expressed genes
                ncgenes=dict42[gene]
                for g1x in ncgenes:
                    g1=g1x.split('==')[0]
                    if g1!=gene:
                        pwys=g1x.split('==')[1].split('$$')
                        if pwy in pwys:
                            v3d[g1]=1
                        else:
                            v4d[g1]=1
            else:
                v3d={}; v4d={}
                            
            #Write to fisherinp if at least 2 genes
            v1=len(list(v1d.keys())); v2=len(list(v2d.keys()))
            v3=len(list(v3d.keys())); v4=len(list(v4d.keys()))
            if v1>=2:
                out2.write('{}\t{}\t{}\t{}\t{}\n'.format(str1,v1,v2,v3,v4))
                out3.write('{}\t{}\t{}\t{}\t{}\n'. \
                           format(str1,list(v1d.keys()),list(v2d.keys()), \
                                  list(v3d.keys()),list(v4d.keys())))
                totcount+=1
            else:
                pass
out1.close(); out2.close(); out3.close()
print ("Total genes in fisherinp: ", totcount)

'''
Do enrichment test COEXP+PWY+   COEXP+PWY-   COEXP-PWY+   COEXP-PWY-
'''

print ("Performing Fisher Exact Test...")
os.system('python2 /home/gdm67/scripts/chengscripts/lsg_Test_Fisher.py {} 1'\
          .format(sys.argv[1]+".coexp.fisherinp"))
    
print ("Making final output file...")
f1=open(sys.argv[1]+".coexp.fisherinp.fisher.pqvalue.sign+pq005.tab", 'r')
#f1=open(sys.argv[1]+".coexp.fisherinp.fisher.pqvalue.sign+pq01.tab", 'r')
line1=f1.readline()
dict7={}
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        sp=tab1[0].split('=='); g1=sp[0]; pwy=sp[1]
        ngenes=tab1[1]; str1='{}|{}'.format(pwy,ngenes)
        #print (sp, g1, pwy, str1)
        #sys.exit()

        if g1 not in dict7:
            dict7[g1]={}
            dict7[g1]['enr']=[str1]
        else:
            if str1 not in dict7[g1]['enr']:
                dict7[g1]['enr'].append(str1)
    line1=f1.readline()
f1.close()

f1=open(sys.argv[2], 'r') #entrez ID file
line1=f1.readline()
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[1]; g2=tab1[2]; pwys=eval(tab1[4])
        if g1 in dlist:

            #Read in alt name
            if g1 not in dict7:
                dict7[g1]={}
                dict7[g1]['alt']=g2
            else:
                if 'alt' not in dict7[g1]:
                    dict7[g1]['alt']=g2

            #Read in PWY info
            if g1 not in dict7:
                dict7[g1]={}
                dict7[g1]['pwy']=pwys
            else:
                if 'pwy' not in dict7[g1]:
                    dict7[g1]['pwy']=pwys
    line1=f1.readline()
f1.close()

out1=open(sys.argv[1]+".coexp.enr.pwys",'w')
out1.write('#python {}\n'.format(' '.join(sys.argv)))
out1.write('#Gene\tAltID\tKnownPWY\tEnrichedPWYS\n')
for gene in dict7:
    altid=dict7[gene]['alt']
    if 'enr' in dict7[gene]:
        enr='; '.join(dict7[gene]['enr'])
    else:
        enr="NA"
    known='; '.join(dict7[gene]['pwy'])
    out1.write('{}\t{}\t{}\t{}\n'.format(gene,altid,known,enr))
out1.close()
print ("See OUT: {}.coexp.enr.pwys".format(sys.argv[1]))
print ("Done!")
