"""-----------------------------------------------------------------------------
  * Script Name: Slyc_geneStructure.py
  * Description: This script will read in a file that contains a list of all 
  BAHD genes and then use that to extract information from a .gff file which
  includes intron and exon count, gene length, and other statistics. 
  * Created By:  Jason David Chobirko
  * Date:        June 28th, 2019
-----------------------------------------------------------------------------"""
# First import all of the necessary modules you need
import sys, os, re

# This list is used to store the open .gff3 files 
fList = []

# This opens and keeps the list of BAHDs needed for the analysis
file = open('Slyc_BAHDs.txt', 'r') 
gList = [line.rstrip("\n") for line in file]

# This creates an output file
newFile = open('slycResults.txt', "w")
newFile.write('gene_ID\tspecies\tfull_gene_length\tcoding_gene_length\texon_length\tintron_length\tintron_count\texon_size\tintron_size\n')

# Idk why the 'r' is necessary, but it made it work...
for f in os.listdir('.'):
	if f.endswith(".gff"):
		fList.append(f)

		
'''
# This iteration will read in each line for each file in fList. Yay! 
for f in fList:
	print("Starting work on " + str(f))
	temp = open(f)
	found = curGene = False
	e1a = e1b = e2a = e2b = gsa = gsb = count = fLen = iLen = eLen = cLen = iNum = 0
	ID = ""
	eList = []; iList = []
	for line in temp:
		if line.startswith("#"):
			if found:
				newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
				count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = [];
			found = False
			continue
		test = re.split('[\t=;\n]', line)
		print("This is test: " + str(test))
		if test[11] in gList and (test[2] == "gene" or test[2] == "mRNA"):
			# If two or more BAHD entries are back-to-back, this should handle that.
			if found:
				newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
				count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = []; 
			ID = test[11]; found = True
			gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb) + 1
			continue
		if found:
			# This is for the other files that do not generously provide introns as a feature. 
			# else:
			if test[2] == "exon":
				# Stand direction is still ignored in these files, so no strand checking needed.
				e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
				num = e1b - e1a + 1; eLen += num
				eList.append(num); count += 1
				# Since introns are between two exons, only once there are at least 2 do we calculate 
				# the necessary information for introns.
				if count >= 2:
					num = (e1a - 1) - (e2b + 1) + 1
					iList.append(num); iLen += num; iNum += 1
			elif test[2] == "CDS":
				# Same as above CDS code
				cLen += int(test[4]) - int(test[3]) + 1
			# 3' and 5' UTRs were triggering this else, which was bad. Fixed?
			elif test[2] == "mRNA" or test[2] == "gene":
				if f == "Azfiliculoides.gene_models.highconfidence_v1.1.gff3" and test[2] == "mRNA":
					continue
				newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
				found = False
		# If we're here it means we hit a gene or mRNA line and it is NOT in the list.
		# Thus we need to print and reset all the variables and carry on our way. 
		else:
			count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; found = False
			continue
newFile.close()
'''			

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
			if found:
				newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
				count = fLen = cLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; #curGene = False
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
					if test[11] in gList:
						# The main method will go in this section which involves intron counting and gene length only.
						ID = test[11]; found = True
						gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb)
					continue
		if test[11] in gList:
			print("Found!")
			# The main method will go in this section which involves intron counting and gene length only.
			ID = test[11]; found = True
			gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb) + 1
		if found:
			if test[2] == "exon":
				# Stand direction is still ignored in these files, so no strand checking needed.
				e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
				num = e1b - e1a + 1; eLen += num
				eList.append(num); count += 1
				# Since introns are between two exons, only once there are at least 2 do we calculate 
				# the necessary information for introns.
				if count >= 2:
					num = (e1a - 1) - (e2b + 1) + 1
					iList.append(num); iLen += num; iNum += 1
			elif test[2] == "CDS":
				# Same as above CDS code
				cLen += int(test[4]) - int(test[3]) + 1			
			'''
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
			'''
		else:
			# We are not within a gene of interest and should reset some variables for safety
			count = fLen = cLen = eLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; curGene = False
			continue
newFile.close()


