"""-----------------------------------------------------------------------------
  * Script Name: Dsal_geneStructure.py
  * Description: This script will read in a file that contains a list of all 
  BAHD genes and then use that to extract information from a .gff file which
  includes intron and exon count, gene length, and other statistics. 
  * Created By:  Jason David Chobirko
  * Date:        July 2nd, 2019
-----------------------------------------------------------------------------"""
# First import all of the necessary modules you need
import sys, os, re

file = open("Dsal.txt", 'r') 
gList = [line.rstrip("\n") for line in file]

newFile = open("Dsal_results.txt", "w")
newFile.write('#gene_ID\tfull_gene_length\tcoding_gene_length\tintron_length\tintron_count\tintron_start\tintron_end\n')

# The code below is much more elegant than code above. Lists are pretty neat from R, so yay.
found = False
e1a = e1b = e2a = e2b = gsa = gsb = count = fLen = iLen = cLen = iNum = 0
ID = ""
iLista = []; iListb = []
for line in open(sys.argv[1]):
	if line.startswith("#"):
		if found:
			newFile.write(ID + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + str(iLista) + '\t' + str(iListb) + '\n')
		found = False
		continue
	test = re.split('[\t=;]', line) 
	if test[9] in gList:
		# The main method will go in this section which involves intron counting and gene length only.
		ID = test[9]; found = True
		gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb)
	if found:			
		if test[2] == "exon":
			e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
			count += 1
			if count >= 2:
				# This is where the distance of the intron will be measured since two are required for an intron's existence
				iLen += min(abs(e1a - e2a), abs(e1a - e2b), abs(e1b - e2a), abs(e1b - e2b)); iNum += 1
				if test[6] == "-":
					iLista.append(e1b - gsa); iListb.append(e2a - gsb)
				else:
					iLista.append(e2b - gsa); iListb.append(e1a - gsb)
		elif test[2] == "CDS":
			cLen += abs(int(test[3]) - int(test[4]))
	else:
		# We are not within a gene of interest and should reset some variables for safety
		count = fLen = cLen = iLen = iNum = 0; iLista = []; iListb = []; ID = ""
		continue

newFile.close()		