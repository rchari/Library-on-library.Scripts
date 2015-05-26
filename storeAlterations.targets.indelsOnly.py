# Author: Raj Chari
# Date: May 21, 2015

# Plasmid pileup: pileup generated from sequencing data of original target plasmid library (used to filter indels from oligo synthesis)
# Samples pileup: mpileup file for all six samples analyzed
# Sample list: comma delimited list of sample names (six in total)
# SP/ST1 -> specify which Cas9 data is being analyzed
# Start -> used to specify the start of the target site (36)
# End -> used to specify the end of the target site (135)
# Alteration summary -> output file summing the different alterations observed (insertions and deletions)
# Reference File -> reference fasta file containing sequence identies of the targets
# Mutation Summary -> output file giving the total # of reads with insertions only, deletions only, both and total read count of target
# Mutation Listing -> list of reconstituted mutated reads broken down per target per sample
# Alteration Count Min -> number of times a particular mutation has to be observed in order to be counted. Could be used to filter spurious indels.


from __future__ import division

import sys
import optparse
import os
import operator
import subprocess
import re
from Bio import SeqIO
from Bio.Seq import Seq

from collections import defaultdict
from collections import Counter
from operator import itemgetter

# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 11:
	raise Exception, "Usage: storeAlterations.targets.indelsOnly.py [plasmidPileup] [samplesPileup] [sampleList] [SP/ST1] [start] [end] [alterationSummary] [referenceFile] [mutationSummary] [mutationListing] [alterationCountMin]"

# store the scoring table
scoringTable = defaultdict(int)
scoreFile = open('Illumina.Scoring.Targets.txt','r')
for line in scoreFile:
	line = line.rstrip('\r\n')
	parts = line.split('\t')
	scoringTable[parts[0]] = int(parts[1])
scoreFile.close()

# complement table
complement = defaultdict(dict)
complement['A'] = 'T'
complement['C'] = 'G'
complement['T'] = 'A'
complement['G'] = 'C'

# function to parse a pileup string
def processPileupString(pileupString,pileupQuality,refBase):
	# set the defaults
	baseList = []
	alterationList = []
	index = 0
	seqIndex = 0
	pileupString = pileupString.upper()
	while index < len(pileupString):
		flagged = 0
		indelStatus = 0
		# case 1: index is a start of a segment
		if pileupString[index]=='^':
			# move two characters to the call
			flagged = 1
			index += 2

		if index < len(pileupString):
			# now check if you have a match or indel
			if flagged==1:
				seqIndex += 1
			else:
				if pileupString[index]=='.':
					baseList.append(refBase)
					seqIndex += 1
				elif pileupString[index]==',':
					baseList.append(complement[refBase])
					seqIndex += 1
				elif pileupString[index]=='*':
					if scoringTable[pileupQuality[seqIndex]] >= 28:
						baseList.append('-')
					else:
						baseList.append(refBase)
					seqIndex += 1
				elif pileupString[index].upper()=='A' or pileupString[index].upper()=='C' or pileupString[index].upper()=='T' or pileupString[index].upper()=='G':
					if scoringTable[pileupQuality[seqIndex]] >= 28:
						baseList.append(pileupString[index])
						alteration = refBase + '>' + pileupString[index].upper() + ':' + str(seqIndex)
						#alterationList.append(alteration)
					else:
						baseList.append(refBase)
					seqIndex += 1
				elif pileupString[index]=='+':
					# check if the next digit is a number
					if pileupString[index+2].isdigit():
						insSize = int(pileupString[index+1:index+3])
						insertion = pileupString[index:index+3+insSize] + ':' + str(seqIndex-1)
						index +=  insSize + 2
					else:
						insSize = int(pileupString[index+1])
						insertion = pileupString[index:index+2+insSize] + ':' + str(seqIndex-1)
						index += insSize + 1


					alterationList.append(insertion)

				elif pileupString[index]=='-':
					# check if the next digit is a number
					if pileupString[index+2].isdigit():
						delStart = index
						delEnd = index + int(pileupString[index+1:index+3]) + 3
						deletion = pileupString[delStart:delEnd] + ':' + str(seqIndex-1)
						index += int(pileupString[index+1:index+3]) + 2
					else:
						delStart = index
						delEnd = index + int(pileupString[index+1]) + 2
						deletion = pileupString[delStart:delEnd] + ':' + str(seqIndex-1)
						index += int(pileupString[index+1]) + 1
					alterationList.append(deletion)
					indelStatus = 1								
		index += 1
	return alterationList

# function to determine if alteration should be counted (i.e. must encompass target site)
def determineAlterationStatus(alteration, position, targetSiteStart, targetSiteEnd):
	countAlteration = 0
	if '+' in alteration or '>' in alteration:
		if position >= targetSiteStart and position <= targetSiteEnd:
			countAlteration = 1
	elif '-' in alteration:
		if alteration[2].isdigit():
			size = int(alteration[1:3])
		else:
			size = int(alteration[1:2])
		if (position >= targetSiteStart and position <= targetSiteEnd) or (position + size >= targetSiteStart and position + size <= targetSiteEnd):
			countAlteration = 1
	return countAlteration

# assigning arguments
plasmidPileup = open(args[0],'r')
samplesPileup = open(args[1],'r')
sampleList = args[2].split(',')
species = args[3]
startIndex = int(args[4])
endIndex = int(args[5])
alterationSummary = open(args[6],'w')
referenceFile = open(args[7])
mutationSummary = open(args[8],'w')
mutationListing = open(args[9],'w')
minCount = int(args[10])

# print statement
print 'Starting file: ' + args[1]

# identify target sites
if species=='SP':
	targetSiteStart = 75
	targetSiteEnd = 97
else:
	targetSiteStart = 73
	targetSiteEnd = 99

# data structures
plasmidLibAlterations = defaultdict(dict)
controlSampleAlterations = defaultdict(dict)
treatedSampleAlterations = defaultdict(dict)
uniqueAlterations = []
depthMatrix = defaultdict(dict)
avgDepthMatrix = defaultdict(dict)
alterationStatusTable = defaultdict(dict)
refDB = defaultdict(str)

# structure to hold the mutated reads -> gene -> sample -> read# -> [alterationList]
mutatedReads = defaultdict(dict)

# store reference sequence
for record in SeqIO.parse(referenceFile,'fasta'):
	refDB[str(record.id)] = str(record.seq)

#print 'Going through plasmid pileup'

# go through plasmid pileup
for line in plasmidPileup:
	line = line.rstrip('\r\n')
	parts = line.split('\t')
	gene = parts[0]
	position = int(parts[1])
	refBase = parts[2].upper()
	# initialize plasmid library
	if gene not in plasmidLibAlterations:
		plasmidLibAlterations[gene] = defaultdict(dict)
	if gene not in controlSampleAlterations:
		controlSampleAlterations[gene] = defaultdict(dict)
	if gene not in treatedSampleAlterations:
		treatedSampleAlterations[gene] = defaultdict(dict)
	if gene not in depthMatrix:
		depthMatrix[gene] = defaultdict(dict)
	if gene not in mutatedReads:
		mutatedReads[gene] = defaultdict(dict)
	if gene not in alterationStatusTable:
		alterationStatusTable[gene] = defaultdict(dict)

	# only grab the positions of interest
	if position >= startIndex and position <= endIndex:

		if parts[1] not in plasmidLibAlterations[gene]:
			plasmidLibAlterations[gene][parts[1]] = defaultdict(int)
		if parts[1] not in controlSampleAlterations[gene]:
			controlSampleAlterations[gene][parts[1]] = defaultdict(int)
		if parts[1] not in treatedSampleAlterations[gene]:
			treatedSampleAlterations[gene][parts[1]] = defaultdict(dict)
		if parts[1] not in depthMatrix[gene]:
			depthMatrix[gene][parts[1]] = defaultdict(dict)
		if parts[1] not in alterationStatusTable[gene]:
			alterationStatusTable[gene][parts[1]] = defaultdict(dict)


		pileupString = parts[4]
		pileupQuality = parts[5]
		alterationSet = processPileupString(pileupString,pileupQuality,refBase)
		# each alteration now holds a position
		alterationList = []
		for alt in alterationSet:
			altParts = alt.split(':')
			alterationList.append(altParts[0])

		myCounter = Counter(alterationList)
		for alteration in myCounter.keys():
			if alteration not in plasmidLibAlterations[gene][parts[1]]:
				plasmidLibAlterations[gene][parts[1]][alteration] = myCounter[alteration]
			else:
				plasmidLibAlterations[gene][parts[1]][alteration] += myCounter[alteration]
			if alteration not in uniqueAlterations:
				uniqueAlterations.append(alteration)
plasmidPileup.close()

#print 'Now going through sample pileup'

# for each sample, you get a matrix gene x position with a listing of the alterations
for line in samplesPileup:
	line = line.rstrip('\r\n')
	parts = line.split('\t')
	gene = parts[0]
	position = int(parts[1])
	refBase = parts[2].upper()	

	# initialize everything
	if gene not in controlSampleAlterations:
		controlSampleAlterations[gene] = defaultdict(dict)
	if gene not in plasmidLibAlterations:
		plasmidLibAlterations[gene] = defaultdict(dict) [position]
	if gene not in treatedSampleAlterations:
		treatedSampleAlterations[gene] = defaultdict(dict)
	if gene not in depthMatrix:
		depthMatrix[gene] = defaultdict(dict)
	if gene not in mutatedReads:
		mutatedReads[gene] = defaultdict(dict)
	if gene not in alterationStatusTable:
		alterationStatusTable[gene] = defaultdict(dict)

	# grab positions of interest
	if position >= startIndex and position <= endIndex:
		if parts[1] not in controlSampleAlterations[gene]:
			controlSampleAlterations[gene][parts[1]] = defaultdict(int)
		if parts[1] not in treatedSampleAlterations[gene]:
			treatedSampleAlterations[gene][parts[1]] = defaultdict(dict)
		if parts[1] not in depthMatrix[gene]:
			depthMatrix[gene][parts[1]] = defaultdict(dict)
		if parts[1] not in plasmidLibAlterations[gene]:
			plasmidLibAlterations[gene][parts[1]] = defaultdict(int)
		if parts[1] not in alterationStatusTable[gene]:
			alterationStatusTable[gene][parts[1]] = defaultdict(dict)

		# only look at situations with all 6 samples
		totalCols = (len(sampleList)*3) + 3
		if len(parts)==totalCols:

			# store alterations for sample 1
			depth = int(parts[3])
			pileupString = parts[4]
			pileupQuality = parts[5]
			alterationSet = processPileupString(pileupString,pileupQuality,refBase)
			# each alteration now holds a position
			alterationList = []
			for alt in alterationSet:
				altParts = alt.split(':')
				alterationList.append(altParts[0])
				if sampleList[0] not in mutatedReads[gene]:
					mutatedReads[gene][sampleList[0]] = defaultdict(list)
				if altParts[1] not in mutatedReads[gene][sampleList[0]]:
					mutatedReads[gene][sampleList[0]][altParts[1]] = []
				# mutation entry
				mutEntry = altParts[0] + ':' + parts[1]
				mutatedReads[gene][sampleList[0]][altParts[1]].append(mutEntry)

			myCounter = Counter(alterationList)
			for alteration in myCounter.keys():
				# identify the alteration span, three cases: insertions, deletions and SNPs
				countAlteration = determineAlterationStatus(alteration, position, targetSiteStart, targetSiteEnd)
				if sampleList[0] not in alterationStatusTable[gene][parts[1]]:
					alterationStatusTable[gene][parts[1]][sampleList[0]] = defaultdict(str)
				if countAlteration==1:
					altParts = alteration.split(':')
					altBaseParts = altParts[0].split('>')
					if '>' not in alteration:
		 				if alteration not in controlSampleAlterations[gene][parts[1]]:
							controlSampleAlterations[gene][parts[1]][alteration] = myCounter[alteration]
						else:
							controlSampleAlterations[gene][parts[1]][alteration] += myCounter[alteration]
						if alteration not in uniqueAlterations:
							uniqueAlterations.append(alteration)


			# adjust the depth
			adjustment = 0
			for c in pileupQuality:
				if scoringTable[c] < 28:
					adjustment += 1

			if sampleList[0] not in avgDepthMatrix:
				avgDepthMatrix[sampleList[0]] = defaultdict(list)

			if gene not in avgDepthMatrix[sampleList[0]]:
				avgDepthMatrix[sampleList[0]][gene] = []


			# add to the depth
			depth = depth - adjustment
			if sampleList[0] not in depthMatrix[gene][parts[1]]:
				depthMatrix[gene][parts[1]][sampleList[0]] = depth
			avgDepthMatrix[sampleList[0]][gene].append(depth)

			# now go through the rest
			x = 1
			while x < len(sampleList):
				# initialize treated
				if sampleList[x] not in treatedSampleAlterations[gene][parts[1]]:
					treatedSampleAlterations[gene][parts[1]][sampleList[x]] = defaultdict(int)

				# avg depth
				if sampleList[x] not in avgDepthMatrix:
					avgDepthMatrix[sampleList[x]] = defaultdict(list)

				if gene not in avgDepthMatrix[sampleList[x]]:
					avgDepthMatrix[sampleList[x]][gene] = []

				# grab the pileup and pileup QL for sample x
				depthPos = ((x + 1) * 3)
				depth = int(parts[depthPos])
				pileupPos = ((x + 1) * 3) + 1
				pqlPos = ((x + 1) * 3) + 2
				alterationSet = processPileupString(parts[pileupPos],parts[pqlPos],refBase)
				# each alteration now holds a position
				alterationList = []
				for alt in alterationSet:
					altParts = alt.split(':')
					alterationList.append(altParts[0])
					if sampleList[x] not in mutatedReads[gene]:
						mutatedReads[gene][sampleList[x]] = defaultdict(list)
					if altParts[1] not in mutatedReads[gene][sampleList[x]]:
						mutatedReads[gene][sampleList[x]][altParts[1]] = []
					# mutation entry
					mutEntry = altParts[0] + ':' + parts[1]
					mutatedReads[gene][sampleList[x]][altParts[1]].append(mutEntry)

				myCounter = Counter(alterationList)
				for alteration in myCounter.keys():
					countAlteration = determineAlterationStatus(alteration, position, targetSiteStart, targetSiteEnd)
					if sampleList[x] not in alterationStatusTable[gene][parts[1]]:
						alterationStatusTable[gene][parts[1]][sampleList[x]] = defaultdict(str)
					if countAlteration==1:
						altParts = alteration.split(':')
						altBaseParts = altParts[0].split('>')
						if '>' not in alteration:						
							if alteration not in treatedSampleAlterations[gene][parts[1]][sampleList[x]]:
								treatedSampleAlterations[gene][parts[1]][sampleList[x]][alteration] = myCounter[alteration]
							else:
								treatedSampleAlterations[gene][parts[1]][sampleList[x]][alteration] += myCounter[alteration]
							if alteration not in uniqueAlterations:
								uniqueAlterations.append(alteration)

				# adjust the depth
				adjustment = 0
				for c in parts[pqlPos]:
					if scoringTable[c] < 28:
						adjustment += 1

				# add to the depth
				depth = depth - adjustment
				if sampleList[x] not in depthMatrix[gene][parts[1]]:
					depthMatrix[gene][parts[1]][sampleList[x]] = depth

				avgDepthMatrix[sampleList[x]][gene].append(depth)

				x += 1

# go through from 36 to 135
header2 = 'Alteration' + '\t' + '\t'.join(sampleList)
alterationSummary.write(header2 + '\n')

# sample specific header
headerPos = 36
headerLine = 'Gene'
allPossiblePositions = []
while headerPos <= 135:
	headerLine = headerLine + '\t' + str(headerPos)
	allPossiblePositions.append(str(headerPos))
	headerPos += 1
alterationCounter = defaultdict(dict)
x = 0
while x < len(sampleList):
	if x==0:
		for alteration in uniqueAlterations:
			count = 0
			# iterate through every gene at every position and look for alteration that matches alteration
			for gene in controlSampleAlterations:
				for position in allPossiblePositions:
					y = 1
					toAdd = 1
					for sampleList[y] in treatedSampleAlterations[gene][position]:
						if alteration in controlSampleAlterations[gene][position] and alteration in treatedSampleAlterations[gene][position][sampleList[y]]:
							toAdd = 0
							break
						y += 1
					if toAdd==1 and alteration in controlSampleAlterations[gene][position] and alteration not in plasmidLibAlterations[gene][position] and controlSampleAlterations[gene][position][alteration] >= minCount:
						count += controlSampleAlterations[gene][position][alteration]
						alterationStatusTable[gene][position][sampleList[0]][alteration] = 'True'
					else:
						alterationStatusTable[gene][position][sampleList[0]][alteration] = 'False'

			if alteration not in alterationCounter:
				alterationCounter[alteration] = defaultdict(int)
			if sampleList[0] not in alterationCounter[alteration]:
				alterationCounter[alteration][sampleList[0]] = count
			else:
				alterationCounter[alteration][sampleList[0]] += count

	else:
		for alteration in uniqueAlterations:
			count = 0
			# iterate through every gene at every position and look for alteration that matches alteration
			for gene in treatedSampleAlterations:
				lineToWrite = gene
				for position in allPossiblePositions:
					if alteration not in controlSampleAlterations[gene][position] and alteration not in plasmidLibAlterations[gene][position] and alteration in treatedSampleAlterations[gene][position][sampleList[x]] and treatedSampleAlterations[gene][position][sampleList[x]][alteration] >= minCount:
						count += treatedSampleAlterations[gene][position][sampleList[x]][alteration]
						alterationStatusTable[gene][position][sampleList[x]][alteration] = 'True'
					else:
						alterationStatusTable[gene][position][sampleList[x]][alteration] = 'False'
			if alteration not in alterationCounter:
				alterationCounter[alteration] = defaultdict(int)
			if sampleList[x] not in alterationCounter[alteration]:
				alterationCounter[alteration][sampleList[x]] = count
			else:
				alterationCounter[alteration][sampleList[x]] += count
				
	x += 1

samplesPileup.close()

#print 'Writing alterations total summary'

for alteration in alterationCounter:
	ltw = alteration
	c = 0
	allZeros = 1
	while c < len(sampleList):
		if alterationCounter[alteration][sampleList[c]] > 0:
			allZeros = 0
		ltw = ltw + '\t' + str(alterationCounter[alteration][sampleList[c]])
		c += 1
	if allZeros==0:
		alterationSummary.write(ltw + '\n')
alterationSummary.close()

# writing the mutation data now
#print 'Writing mutation data'
x = 0
header = 'Gene'
while x < len(sampleList):
	header = header + '\t' + sampleList[x] + '.InsReads\t' + sampleList[x] + '.DelReads\t' + sampleList[x] + '.MixedReads\t' + sampleList[x] + '.TotalReads'
	x += 1
mutationSummary.write(header + '\n')

# go through the table
for gene in mutatedReads:

	# write the gene header
	mutationListing.write('#####GENE#####\n' + gene + '\n')

	# write the reference sequence
	mutationListing.write('>ReferenceSequence\n' + refDB[gene] + '\n')

	# tally the mutations
	sampleNumber = 0
	summToWrite = gene
	while sampleNumber < len(sampleList):
		# grab the max read count (depth)
		if len(avgDepthMatrix[sampleList[sampleNumber]][gene]) > 0:
			totalReads = max(avgDepthMatrix[sampleList[sampleNumber]][gene])
		else:
			totalReads = 0

		# initialize counts
		insCount = 0
		delCount = 0
		snpCount = 0
		mixedCount = 0

		# write the sample name
		mutationListing.write('Sample name: ' + sampleList[sampleNumber] + '\n')

		if sampleList[sampleNumber] in mutatedReads[gene]:
			for readID in mutatedReads[gene][sampleList[sampleNumber]]:
				clearedAlterations = []
				for alteration in mutatedReads[gene][sampleList[sampleNumber]][readID]:
					altParts = alteration.split(':')
					#print 'Alteration: ' + alteration + ', Status: ' + alterationStatusTable[gene][altParts[1]][sampleList[sampleNumber]][altParts[0]]
					if alterationStatusTable[gene][altParts[1]][sampleList[sampleNumber]][altParts[0]]=='True':
						clearedAlterations.append(alteration)

				#print 'Size of cleared alterations: ' + str(len(clearedAlterations))

				# now update the counts
				if len(clearedAlterations) > 1:
					mixedCount += 1					
				elif len(clearedAlterations)==1: 
					if '+' in clearedAlterations[0]:
						insCount += 1
					elif '-' in clearedAlterations[0]:
						delCount += 1
					elif '>' in clearedAlterations[0]:
						snpCount += 1


				if len(clearedAlterations) > 0:
					refStart = 0
					mutatedSequence = ''
					for alteration in clearedAlterations:
						# alteration : position
						altParts = alteration.split(':')
						altPosition = int(altParts[1])

						# if it's a SNP, it's the same right to pos - 1
						if '>' in altParts[0]:
							bases = altParts[0].split('>')
							mutatedSequence = mutatedSequence + refDB[gene][refStart:altPosition-1] + bases[1]

							# update refStart
							refStart = altPosition
						else:
							# get the deletion size
							event = altParts[0]
							if event[2].isdigit():
								eventSize = int(event[1:3])
								eventSeq = event[3:]
							else:
								eventSize = int(event[1:2])
								eventSeq = event[2:]

							if '-' in event:
								mutatedSequence = mutatedSequence + refDB[gene][refStart:altPosition] + eventSize * '-'
								# update refStart
								refStart = altPosition + eventSize
							else:
								mutatedSequence = mutatedSequence + refDB[gene][refStart:altPosition] + eventSeq
								# update refStart
								refStart = altPosition
					# add the remaining sequence
					if refStart < 170:
						mutatedSequence = mutatedSequence + refDB[gene][refStart:]

					# write mutated sequence accordingly
					mutationListing.write('>' + str(readID) + ':' + ';'.join(clearedAlterations) + '\n' + mutatedSequence + '\n')

		# if the count is -1 (no reads, set to 0)
		if totalReads==-1:
			totalReads = 0

		summToWrite = summToWrite + '\t' + str(insCount) + '\t' + str(delCount) + '\t' + str(mixedCount) + '\t' + str(totalReads)
		sampleNumber += 1

	# write to summary file
	mutationSummary.write(summToWrite + '\n')

mutationSummary.close()
mutationListing.close()

