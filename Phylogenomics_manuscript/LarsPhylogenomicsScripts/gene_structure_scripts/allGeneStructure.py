"""-----------------------------------------------------------------------------
  * Script Name: allGeneStructure.py
  * Description: This script will read in a file that contains a list of all 
  BAHD genes and then use that to extract information from a .gff file which
  includes intron and exon count, gene length, and other statistics. 
  * Created By:  Jason David Chobirko
  * Date:        July 5th, 2019
-----------------------------------------------------------------------------"""
# First import all of the necessary modules you need
import sys, os, re
#from itertools import zip_longest

# This list is used to store the open .gff3 files 
fList = []

# This opens and keeps the list of BAHDs needed for the analysis
file = open('allBAHDs.txt', 'r') 
gList = [line.rstrip("\n") for line in file]

# This creates an output file
newFile = open('allResults.txt', "w")
newFile.write('gene_ID\tspecies\tfull_gene_length\tcoding_gene_length\texon_length\tintron_length\tintron_count\texon_size\tintron_size\n')

# This is for convience so I don't have to move stuff anywhere. Also why is the 'r' needed?
os.chdir(r"C:\Users\Jason Chobirko\Desktop\Phytozome_GFFs")

# Idk why the 'r' is necessary, but it made it work...
for f in os.listdir('.'):
	if f.endswith(".gff3"):
		fList.append(f)

# This iteration will read in each line for each file in fList. Yay! 
for f in fList:
	print("Starting work on " + str(f))
	temp = open(f)
	found = curGene = False
	e1a = e1b = e2a = e2b = gsa = gsb = count = fLen = iLen = cLen = eLen = iNum = 0
	ID = ""
	eList = []; iList = []
	for line in temp:
		# The code below is much more elegant than code above. Lists are pretty neat from R, so yay.	
		if line.startswith("#"):
			found = False
			continue
		test = re.split('[\t=;]', line)
		if test[2] == "gene":
			#print("gene was found at:" + str(test))
			curGene = True
			if found:
				newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
			found = False
			count = fLen = cLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; #curGene = False
		if test[2] == "mRNA":
			if curGene:
				curGene = False
			else:
				if found:
					newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
					found = False; 
					count = fLen = cLen = eLen = iLen = iNum = 0; eList = []; iList = [];
					if test[9] in gList:
						# The main method will go in this section which involves intron counting and gene length only.
						ID = test[9]; found = True
						gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb)
					continue
		if test[9] in gList:
			# The main method will go in this section which involves intron counting and gene length only.
			ID = test[9]; found = True
			gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb) + 1
		if found:
			if test[2] == "exon":
				e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
				num = e1b - e1a + 1; eLen += num
				eList.append(num); count += 1
				if count >= 2:
					# This is where the distance of the intron will be measured since two are required for an intron's existence
					num = 0; iNum += 1
					if test[6] == "-":
						num = (e2a - 1) - (e1b + 1) + 1
						iList.append(num)
						iLen += num
					else:
						num = (e1a - 1) - (e2b + 1) + 1
						iList.append(num)
						iLen += num
			elif test[2] == "CDS":
				cLen += abs(int(test[3]) - int(test[4])) + 1
		else:
			# We are not within a gene of interest and should reset some variables for safety
			count = fLen = cLen = eLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; curGene = False
			continue
newFile.close()