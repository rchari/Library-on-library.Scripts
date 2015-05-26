# Author: Raj Chari
# Date: May 21st, 2015
# script to take a pileup and split it into one file for each gene/target

from __future__ import division

import sys
import optparse
import os
import operator
import subprocess

from collections import defaultdict
from operator import itemgetter

# parse the command line
parser = optparse.OptionParser()
(options,args) = parser.parse_args()
if len(args) != 2:
	raise Exception, "Usage: splitPileup.py [mPileupFile] [sampleNames]\n"
# assign file handles
infile = open(args[0],'r')
s1 = args[1]

# go through the pileup file and create a separate pilueup file for every gene
currSample = ''
lw = 0
for line in infile:
	line = line.rstrip('\r\n')
	parts = line.split('\t')
	if int(parts[1]) >= 36 and int(parts[1]) <= 135:	
		if currSample=='':
			filename = s1 + '.' + parts[0] + '.' + '36.135.Pileup.txt'
			outfile = open(filename,'w')

		elif parts[0] != currSample:
			#print 'Done sample: ' + currSample + ', Wrote this many lines: ' + str(lw)
			outfile.close()
			# check here if 135 lines were written, if not, delete the file
			if lw!=100:
				os.remove(filename)
			lw = 0
			filename = s1 + '.' + parts[0] + '.' + '36.135.Pileup.txt'
			outfile = open(filename,'w')

		outfile.write(line + '\n')
		currSample = parts[0]
		lw += 1
outfile.close()
infile.close()
