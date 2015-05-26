# Author: Raj Chari
# Date: May 21st, 2015

# script to extract separate pileups for each gene

from __future__ import division

import sys
import optparse
import os
import operator
import re

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 2:
    raise Exception, "python extractFromPileupToGene.py [pileupFile] [geneFile]"

# file handles
pileup = open(args[0],'r')
regions = open(args[1],'r')

# group filename
gfn = args[0].replace('.mPileup.txt','')

# go through the gene file
regionDict = defaultdict(list)
pileupDict = defaultdict(list)
regionCount = 0

for line in regions:
	line = line.rstrip('\r\n')
	parts = line.split('\t')
	regionDict[parts[0]].append(line)
	regionCount += 1
regions.close()


# go through pileup
currentRegion = ''
lineCount = 0
for line in pileup:
	line = line.rstrip('\r\n')
	parts = line.split('\t')

	if lineCount % 1000==0:
		print 'Line count: ' + str(lineCount) + ', Chromosome: ' + parts[0]

	# find the right chromsome
	ROIs = regionDict[parts[0]]

	# check if the current region is still good
	if currentRegion != '':
		currParts = currentRegion.split('\t')
		if int(parts[1]) >= int(currParts[1]) and int(parts[1]) <= int(currParts[2]) and currParts[0]==parts[0]:
			output.write(line + '\n')
		else:
			if currentRegion in regionDict[parts[0]]:
				regionDict[parts[0]].remove(currentRegion)
				regionCount = regionCount - 1
				print 'Regions left: ' + str(regionCount)
				output.close()
			currentRegion = ''
			if len(ROIs) > 0:
				for region in ROIs:
					regionParts = region.split('\t')
					if int(parts[1]) >= int(regionParts[1]) and int(parts[1]) <= int(regionParts[2]):
						fileName = gfn + '.' + regionParts[3] + '.txt'
						output = open(fileName,'w')
						output.write(line + '\n')
						currentRegion = region
						break
	else:
		if len(ROIs) > 0:
			for region in ROIs:
				regionParts = region.split('\t')
				if int(parts[1]) >= int(regionParts[1]) and int(parts[1]) <= int(regionParts[2]):
					fileName = gfn + '.' + regionParts[3] + '.txt'
					output = open(fileName,'w')
					output.write(line + '\n')
					currentRegion = region
					break

	if regionCount==0:
		break

	lineCount += 1

pileup.close()

# go through pileup dictionary to write files
for gene in pileupDict:

	for line in pileupDict[gene]:
		output.write(line + '\n')
	output.close()
