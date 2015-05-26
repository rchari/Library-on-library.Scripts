# Author: Raj Chari
# Date: May 21st, 2015

# script to merge the mutation listings, mutation summaries, and the alteration stats
# input
# after running storeAlterations, each target will have an alteration summary, mutation summary and mutation listing
# create a file that lists all the files for each type of output.  In total, there will be list files that will be used as input
# sample list is comma delimited (i.e. WT1,WT2,WT3,WT4,WT5,WT6)
# group header will mark every output file created with <group header>.*

from __future__ import division

import sys
import optparse
import os
import operator
import subprocess
import re
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from collections import defaultdict
from collections import Counter
from operator import itemgetter

# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 5:
	raise Exception, "Usage: mergeFiles.py [List.MutationSummary] [List.MutationListing] [List.AlterationStats] [sampleList] [groupHeader]\n"

# assign the file handles
msList = open(args[0],'r')
mlList = open(args[1],'r')
asList = open(args[2],'r')
sampleList = args[3]
groupName = args[4]

# use the group header to determine the header names
insCol = '.InsReads'
delCol = '.DelReads'
mixCol = '.MixedReads'
totalCol = '.TotalReads'
colNameList = []
mutSummaryDB = defaultdict(dict)

sampleNameList = sampleList.split(',')

x = 0
while x < len(sampleNameList):
	colNameI = sampleNameList[x] + insCol 
	colNameD = sampleNameList[x] + delCol
	colNameM = sampleNameList[x] + mixCol
	colNameT = sampleNameList[x] + totalCol
	colNameList.append(colNameI)
	colNameList.append(colNameD)
	colNameList.append(colNameM)
	colNameList.append(colNameT)	
	x += 1

# debug
#print colNameList

# go through each list; first the msList
fileNumber = 0
for myFile in msList:
	myFile = myFile.rstrip('\r\n')
	lc = 0
	infile = open(myFile,'r')

	# get the ending of the file name
	nameParts = myFile.split('.')
	ending = '.'.join(nameParts[-4:])

	for line in infile:
		line = line.rstrip('\r\n')
		if lc==0:
			colOrder = line
		else:
			colNames = colOrder.split('\t')
			data = line.split('\t')
			mutSummaryDB[data[0]] = defaultdict(int)
			c = 1
			while c < len(colNames):
				mutSummaryDB[data[0]][colNames[c]] = int(data[c])
				#print 'Col Name: ' + colNames[c] + ', Value: ' + data[c]
				c += 1
		lc += 1
	infile.close()
	fileNumber += 1

# file handles for output
msFile = groupName + '.MutationSummary.' + ending
mlFile = groupName + '.MutationListing.' + ending
asFile = groupName + '.AlterationSummary.' + ending
msf = open(msFile,'w')
mlf = open(mlFile,'w')
asf = open(asFile,'w')

# write the mutation summary
headerLine = '\t'.join(colNameList)
msf.write('Gene' + '\t' + headerLine + '\n')
for gene in mutSummaryDB:
	c = 0
	ltw = gene
	while c < len(colNameList):
		ltw = ltw + '\t' + str(mutSummaryDB[gene][colNameList[c]])
		c += 1
	msf.write(ltw + '\n')

# next the ml list
for myFile in mlList:
	myFile = myFile.rstrip('\r\n')
	infile = open(myFile,'r')
	for line in infile:
		mlf.write(line)
	infile.close()

# finally write the asFile
asData = defaultdict(dict)
fileNumber = 0
for myFile in asList:
	myFile = myFile.rstrip('\r\n')
	lc = 0
	infile = open(myFile,'r')
	for line in infile:
		line = line.rstrip('\r\n')
		if lc==0:
			colNames = line.split('\t')
		else:
			parts = line.split('\t')
			if parts[0] not in asData:
				asData[parts[0]] = defaultdict(int)
			c = 1
			while c < len(parts):
				if colNames[c] not in asData[parts[0]]:
					asData[parts[0]][colNames[c]] = int(parts[c])
				else:
					asData[parts[0]][colNames[c]] += int(parts[c])
				c += 1
		lc += 1
	infile.close()

# write alterations
headerLine = '\t'.join(colNames)
asf.write(headerLine + '\n')
for alteration in asData:
	ltw = alteration
	c = 1
	while c < len(colNames):
		ltw = ltw + '\t' + str(asData[alteration][colNames[c]])
		c += 1
	asf.write(ltw + '\n')

# close all files
msList.close()
mlList.close()
asList.close()
msf.close()
mlf.close()
asf.close()
