"""-----------------------------------------------------------------------------
  * Script Name: algaeGeneStructure.py
  * Description: This is a version of the general "allGeneStructure.py" script
  that is dealing with .gff(3) files from the data-gathering of algae. Species 
  include green, red, brown and other various algae. The files are kinda weird
  so this script intends to handle them all in one simple stroke. Strands are 
  followed, so that is needed to keep in mind. 
  * Created By:  Jason David Chobirko
  * Date:        July 12th, 2019
-----------------------------------------------------------------------------"""
# First import all of the necessary modules you need
import sys, os, re
#from itertools import zip_longest

# This list is used to store the open .gff3 files 
fList = []

# This opens and keeps the list of BAHDs needed for the analysis
file = open('algaeBAHDs.txt', 'r') 
gList = [line.rstrip("\n") for line in file]

# This creates an output file
newFile = open('algaeResults.txt', "w")
newFile.write('gene_ID\tspecies\tfull_gene_length\tcoding_gene_length\texon_length\tintron_length\tintron_count\texon_size\tintron_size\n')

# This is for convience so I don't have to move stuff anywhere. Also why is the 'r' needed?
os.chdir(r"C:\Users\Jason Chobirko\Desktop\Phytozome_GFFs\Algae\Test")

# Idk why the 'r' is necessary, but it made it work...
for f in os.listdir('.'):
	if f.endswith(".gff3") or f.endswith(".gff"):
		fList.append(f)

# This iteration will read in each line for each file in fList. Yay! 
for f in fList:
	print("Starting work on " + str(f))
	#"""
	temp = open(f)
	found = curGene = False
	e1a = e1b = e2a = e2b = gsa = gsb = count = fLen = iLen = eLen = cLen = iNum = sta = stb = ena = enb = 0
	ID = ""
	eList = []; iList = []
	for line in temp:
		if line.startswith("#"):
			if found:
				newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
				count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = [];
			found = False
			continue
			
		test = re.split('[\t=;\n"]', line)
		
		#print(str(test))

		### This section is to deal with the specific location of various .gff3 files
		if f == "Cpapaya_113_ASGPBv0.4.gene_exons.gff3":
			if test[11] in gList and (test[2] == "gene" or test[2] == "mRNA"):
				#print("This BAHD was found: " + test[11])
				#print("This is the current eLen: " + str(eLen) + str(found))
				#print("This is the NEW eLen: " + str(eLen))
				ID = test[11]; found = True
				gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb) + 1
				continue
			if found:
				if test[2] == "exon":
					# Stand direction is NOT ignored here, so strand checking is needed.
					e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
					num = e1b - e1a + 1; eLen += num
					eList.append(num); count += 1
					# Since introns are between two exons, only once there are at least 2 do we calculate 
					# the necessary information for introns.
					if count >= 2:
						if test[6] == "-":
							num = (e2a - 1) - (e1b + 1) + 1
							iList.append(num); iLen += num; iNum += 1
						else:
							num = (e1a - 1) - (e2b + 1) + 1
							iList.append(num); iLen += num; iNum += 1
				elif test[2] == "CDS":
					# Same as above CDS code
					cLen += int(test[4]) - int(test[3]) + 1
				# 3' and 5' UTRs were triggering this else, which was bad. Fixed?
				elif test[2] == "mRNA" or test[2] == "gene":
					#print("This was hit for ID: " + str(test[9]))
					newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
					count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; found = False
					a = 0
			# If we're here it means we hit a gene or mRNA line and it is NOT in the list.
			# Thus we need to print and reset all the variables and carry on our way. 
			else:
				#newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
				count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; found = False
				continue
			
		# This is for the majority of 'non-weird' algae files 
		else:
		
			if f == "Mcommoda.gff3" or f == "Ota2ri.gff3" or f == "Picochlorum.gff3":
				if (f[:4] + "_" + test[9]) in gList and (test[2] == "gene" or test[2] == "mRNA"):
					found = True; ID = test[9]; gsa = int(test[3]); gsb = int(test[4]); fLen += abs(gsa - gsb) + 1
					continue
				if found:
					if test[2] == "exon":
						# Stand direction is NOT ignored here, so strand checking is needed.
						e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
						num = e1b - e1a + 1; eLen += num
						eList.append(num); count += 1
						# Since introns are between two exons, only once there are at least 2 do we calculate 
						# the necessary information for introns.
						if count >= 2:
							if test[6] == "-":
								num = (e2a - 1) - (e1b + 1) + 1
								iList.append(num); iLen += num; iNum += 1
							else:
								num = (e1a - 1) - (e2b + 1) + 1
								iList.append(num); iLen += num; iNum += 1
					elif test[2] == "CDS":
						# Same as above CDS code
						cLen += int(test[4]) - int(test[3]) + 1
					elif test[2] == "gene": # test[2] == "mRNA" or
						newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
						found = False
				else:
					count = fLen = eLen = cLen = iLen = iNum = 0; eList = []; iList = []; ID = ""; found = False
					continue
			# The files below have their ID match in EVERY feature, which needs to be handled differently
			else: 
				if (f[:4] + "_" + test[9]) in gList:
					found = True; ID = test[9]
					if test[2] == "exon":
						e2a = e1a; e1a = int(test[3]); e2b = e1b; e1b = int(test[4])
						num = e1b - e1a + 1; eLen += num
						eList.append(num); count += 1
						# Since introns are between two exons, only once there are at least 2 do we calculate 
						# the necessary information for introns.
						if count >= 2:
							if test[6] == "-":
								num = (e2a - 1) - (e1b + 1) + 1
								iList.append(num); iLen += num; iNum += 1
							else:
								num = (e1a - 1) - (e2b + 1) + 1
								iList.append(num); iLen += num; iNum += 1
					elif test[2] == "CDS":
						# Same as above CDS code
						cLen += int(test[4]) - int(test[3]) + 1
					elif test[2] == "start_codon":
						sta = int(test[3]); stb = int(test[4])
					elif test[2] == "stop_codon":
						ena = int(test[3]); enb = int(test[4])
				else:
					if found:
						fLen = abs(sta - enb) + 1
						newFile.write(ID + '\t' + "{0}".format(f)[:4] + '\t' + str(fLen) + '\t' + str(cLen) + '\t' + str(eLen) + '\t' + str(iLen) + '\t' + str(iNum) + '\t' + ','.join(map(str, eList)) + '\t' + ','.join(map(str, iList)) + '\n')
					count = fLen = eLen = cLen = iLen = iNum = sta = stb = ena = enb = 0; eList = []; iList = []; ID = ""; found = False
					continue
newFile.close()