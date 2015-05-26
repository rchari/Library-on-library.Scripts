# Scripts used in the study by Chari et al, Nature Methods, 2015

# Description of each file:

storeAlterations.targets.indelsOnly.py - python script that outputs mutation rates from the integrated target dataset.  Pileup files should be split on a per target basis to speed up run time. Prior to running this script, run the "splitPileup.py" on the aggregate pileup file containing all six samples.

splitPileup.py - python script to split a pileup file into individual files for each target site

Illumina.Scoring.Targets.txt - score file used by the storeAlterations script to use base quality scores for filtering





