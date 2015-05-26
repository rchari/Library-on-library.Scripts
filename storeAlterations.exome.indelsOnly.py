# Author: Raj Chari
# Date: May 21st, 2015


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
if len(args) != 3:
	raise Exception, "Usage: storeAlterations.exome.indelsOnly.py [pileupList] [sampleList] [mutationSummary]"

# complement table
complement = defaultdict(dict)
complement['A'] = 'T'
complement['C'] = 'G'
complement['T'] = 'A'
complement['G'] = 'C'

# assigning arguments
samplesPileup = open(args[0],'r')
sampleList = args[1].split(',')
mutationSummary = open(args[2],'w')

# header lines
mutationSummary.write('GeneName\tRefSeq\tControl.Mutations\tControl.Total\tTreated.Mutations\tTreated.Total\n')

# function to parse a pileup string
def processPileupString(pileupString,pileupQuality,refBase):
	# set the defaults
	baseList = []
	alterationList = []
	index = 0
	seqIndex = 0
	pileupString = pileupString.upper()
	while index < len(pileupString):

		#print 'Char: ' + pileupString[index] + ', Size of baselist: ' + str(len(baseList)) + ', SeqIndex: ' + str(seqIndex)

		# case 1: index is a start of a segment
		if pileupString[index]=='^':
			# move two characters to the call
			index += 2

		if index < len(pileupString):

			if pileupString[index]=='.':
				baseList.append(refBase)
				seqIndex += 1
			elif pileupString[index]==',':
				baseList.append(complement[refBase])
				seqIndex += 1
			elif pileupString[index]=='*':
				baseList.append('D')
				seqIndex += 1
			elif pileupString[index].upper()=='A' or pileupString[index].upper()=='C' or pileupString[index].upper()=='T' or pileupString[index].upper()=='G':
				if ord(pileupQuality[seqIndex]) - 64 >= 28:
					alteration = refBase + '>' + pileupString[index].upper()
					baseList.append(alteration)
				else:
					baseList.append(refBase)
				seqIndex += 1
			elif pileupString[index]=='+':
				# check if the next digit is a number
				if pileupString[index+2].isdigit():
					insSize = int(pileupString[index+1:index+3])
					insertion = pileupString[index:index+3+insSize]
					index +=  insSize + 2
				else:
					insSize = int(pileupString[index+1])
					insertion = pileupString[index:index+2+insSize]
					index += insSize + 1
				baseList[seqIndex-1] = baseList[seqIndex-1] + insertion

			elif pileupString[index]=='-':
				# check if the next digit is a number
				if pileupString[index+2].isdigit():
					delStart = index
					delEnd = index + int(pileupString[index+1:index+3]) + 3
					deletion = pileupString[delStart:delEnd]
					index += int(pileupString[index+1:index+3]) + 2
				else:
					delStart = index
					delEnd = index + int(pileupString[index+1]) + 2
					deletion = pileupString[delStart:delEnd]
					index += int(pileupString[index+1]) + 1
				baseList[seqIndex-1] = baseList[seqIndex-1] + deletion							
		index += 1
	return baseList

# function to determine relevant alteration
def determineAlterationStatus(alteration, position, targetSiteStart, targetSiteEnd):
	countAlteration = 0
	if '+' in alteration or '>' in alteration:
		if position >= targetSiteStart and position <= targetSiteEnd:
			countAlteration = 1
	elif '-' in alteration:
		if alteration[3].isdigit():
			size = int(alteration[2:4])
		else:
			size = int(alteration[2:3])
		if (position >= targetSiteStart and position <= targetSiteEnd) or (position + size >= targetSiteStart and position + size <= targetSiteEnd):
			countAlteration = 1
	return countAlteration


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

print 'Now going through sample pileup'

# for each sample, you get a matrix gene x position with a listing of the alteration
for infile in samplesPileup:
	infile = infile.rstrip('\r\n')
	infileHandle = open(infile,'r')

	offSetS1 = 0
	offSetS2 = 0
	seqListS1 = []
	seqListS2 = []
	gene = infile
	minPosition = 99999999999999999
	refSeq = ''
	altListS1 = []
	altListS2 = []

	print 'Starting file: ' + infile

	for line in infileHandle:

		# break each part of the pileup
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		chromosome = parts[0]
		position = int(parts[1])
		refBase = parts[2].upper()
		refSeq = refSeq + refBase

		# need this for relative position
		if position < minPosition:
			minPosition = position

		# break down base lists
		baseListS1 = processPileupString(parts[4],parts[5],refBase)
		#print 'Position: ' + str(position)
		baseListS2 = processPileupString(parts[7],parts[8],refBase)

		# append to the list
		x = 0
		while x < len(baseListS1):
			alteration = baseListS1[x] + ':' + str(position) + ';'
			if x + offSetS1 + 1 > len(seqListS1):		
				seqListS1.append(alteration)
			else:
				seqListS1[x+offSetS1] = seqListS1[x+offSetS1] + alteration

			# add to altList
			if '+' in baseListS1[x] or '-' in baseListS1[x]:
				if alteration not in altListS1:
					altListS1.append(alteration)

			x += 1
		y = 0
		while y < len(baseListS2):
			alteration = baseListS2[y] + ':' + str(position) + ';'
			if y + offSetS2 + 1 > len(seqListS2):		
				seqListS2.append(alteration)
			else:
				seqListS2[y+offSetS2] = seqListS2[y+offSetS2] + alteration

			# add to altList
			if '+' in baseListS2[y] or '-' in baseListS2[y]:
				if alteration not in altListS2:
					altListS2.append(alteration)
			y += 1

		# calculate the offsets
		offSetS1 = parts[4].count('$')
		offSetS2 = parts[7].count('$')

	# define target start and target end
	if 'ST1.' in infile:
		targetSiteStart = minPosition + 100
		targetSiteEnd = targetSiteStart + 27
	else:
		targetSiteStart = minPosition + 100
		targetSiteEnd = targetSiteStart + 23


	# go through sample 1:
	s1MutCount = 0
	s1TotalCount = 0
	s2MutCount = 0
	s2TotalCount = 0

	for sequence in seqListS1:
		elements = sequence.split(';')
		spanTargetSite = defaultdict(str)
		tss = targetSiteStart
		tse = targetSiteEnd
		while tss <= tse:
			spanTargetSite[str(tss)] = 'N'
			tss += 1

		for element in elements:
			if '+' in element or '-' in element:
				elemParts = element.split(':')
				countAlteration = determineAlterationStatus(elemParts[0], int(elemParts[1]), targetSiteStart, targetSiteEnd)
				element = element + ';'
				if element not in altListS2 and countAlteration==1:
				#print 'Control: ' + element
					s1MutCount += 1
			# count read if it spans target site
			if element != '':
				elemParts = element.split(':')
				spanTargetSite[elemParts[1]] = 'Y'

		countForTotal = 1
		for position in spanTargetSite:
			if spanTargetSite[position]=='N':
				countForTotal = 0
		if countForTotal==1:
			s1TotalCount += 1

	for sequence in seqListS2:
		elements = sequence.split(';')
		spanTargetSite = defaultdict(str)
		tss = targetSiteStart
		tse = targetSiteEnd
		while tss <= tse:
			spanTargetSite[str(tss)] = 'N'
			tss += 1
		for element in elements:
			if '+' in element or '-' in element:
				elemParts = element.split(':')
				countAlteration = determineAlterationStatus(elemParts[0], int(elemParts[1]), targetSiteStart, targetSiteEnd)
				element = element + ';'
 				if element not in altListS1 and countAlteration==1:
				#print 'Treated: ' + element
					s2MutCount += 1
			# count read if it spans target site
			if element != '':
				elemParts = element.split(':')
				spanTargetSite[elemParts[1]] = 'Y'
		
		countForTotal = 1
		for position in spanTargetSite:
			if spanTargetSite[position]=='N':
				countForTotal = 0
		if countForTotal==1:
			s2TotalCount += 1

	mutationSummary.write(infile + '\t' + refSeq + '\t' + str(s1MutCount) + '\t' + str(s1TotalCount) + '\t' + str(s2MutCount) + '\t' + str(s2TotalCount) + '\n')

mutationSummary.close()
