import sys

print ("INP1: Orthofinder groups file (Orthogroups.txt)")
print ("INP2: List of files where speciesName is first part before _ e.g. Acoe_bahd.fa")
print ("INP3: BLASTP top match")


def read_orthogroups(fn, D1, D2):    
    file1=open(fn,'r')
    line1=file1.readline()    
    while line1:
        if line1.startswith('#'):
            pass
        else:
            if line1.startswith('OG'):
                sp=line1.strip().split(); genes=sp[1:]
                og=sp[0].split(':')[0]                
                for gene in genes:
                    if gene not in D1:
                        D1[gene]=og
                    else:
                        print ("Gene repeat: ", gene)
                        sys.exit()
                D2[og]=genes
        line1=file1.readline()
    file1.close()
    return D1,D2

def read_names(flist, D):
    for line in flist:
        sp=line.strip().split('\t')
        species=sp[0].split('_')[0]     
        sprange=sp[2]
        D[species]=sprange
    return D

#Read orthogroups file
dict1={}; dict2={}
dict1,dict2=read_orthogroups(sys.argv[1], dict1, dict2)

#Read species names files and get species prefix
filelist=open(sys.argv[2],'r').readlines()
sdict={}
sdict=read_names(filelist, sdict)

#Process Results
file1=open(sys.argv[3],'r') #BLASTP topmatch
out1=open(sys.argv[3]+".og",'w')
out2=open(sys.argv[3]+".og.sp",'w')
line1=file1.readline()
dict3={}; x=0; y=0; z=0
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[0]; g2=tab1[1]; idt=tab1[2]; mlen=tab1[3]
        dict3[g1]=1
        if g2 in dict1:
            og=dict1[g2]
            genes=dict2[og]
            out1.write('{}\t{}\t{}\t{}\n'.format(g1,g2,og,','.join(genes)))
            y+=1

            #Get species ranges
            dict4={}
            for gene in genes:
                abbr=gene.split('_')[0]
                dict4[abbr]=1
                
            newlist=[]; rlist=[]
            for species in sdict:                
                if species in dict4:
                    newlist.append(species)
                    srange=sdict[species]
                    rlist.append(srange)

            sline=','.join(newlist)
            rline=','.join(rlist)
            out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(g1,g2,og,idt,mlen,sline,rline))
                
        else:
            out1.write('{}\t{}\tNA\tNA\n'.format(g1,g2))
            print (g1, g2)
            z+=1
    line1=file1.readline()
file1.close(); out1.close(); out2.close()
print ("# of genes: ", len(dict3.keys()))
print ("# of genes with OG: ", y)
print ("# of genes without OG: ", z)


adict={'red_algae':1,'green_algae':2,'charophytic_green_algae':3, \
       'liverworts':4,'mosses':5,'hornworts':6,'lycophytee':7,  'ferns':8, \
       'gymnosperms':9,'angiosperms':10,'dicots':11,'monocots':12}

bdict={1:'red_algae',2:'green_algae',3:'charophytic_green_algae', \
       4:'liverworts',5:'mosses',6:'hornworts',7:'lycophytee',8:'ferns', \
       9:'gymnosperms',10:'angiosperms',11:'dicots',12:'monocots'}

file1=open(sys.argv[3]+".og.sp",'r')
out1=open(sys.argv[3]+".og.sp.range",'w')
line1=file1.readline()
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        try:
            trange=tab1[6].split(',')
            #Get dicot-monocot out first
            taxs=list(set(trange))
            if len(taxs)==2 and 'dicots' in taxs and 'monocots' in taxs:
                finrange='dicot-monocot'
            else:        
                vlist=[]
                for tax in trange:
                    val=adict[tax]
                    if vlist==[]:
                        vlist.append(val)
                    else:
                        if val<vlist[0]:
                            vlist[0]=val
                finval=vlist[0]
                finrange=bdict[finval]
        except:
            finrange='Taxus'
               
        out1.write('{}\t{}\n'.format(line1.strip(),finrange))
    line1=file1.readline()
file1.close(); out1.close()
print ("Done!")


            
                
