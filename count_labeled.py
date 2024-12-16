#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

"""
Count number of labeled PSMs, peptides, and proteins

Arguments:
-psm: psm.txt file
-pep: pep.txt file
-pro: pro.txt file
-l: .LabelPCTcount.txt file 
-th: enrichment threshold to use 

Output: 
Sipros labeled search summary stats
- Total number of PSMs
- Number of labeled PSMs
- Total number of peptides
- Number of labeled peptides
- Total number of proteins
- Number of labeled proteins
"""

parser = argparse.ArgumentParser()
parser.add_argument("-psm", "--psmTxtFile")
parser.add_argument("-pep", "--pepTxtFile")
parser.add_argument("-pro", "--proTxtFile")
parser.add_argument("-l", "--lblPCT")
parser.add_argument("-th", "--threshold")
args = parser.parse_args()


### getSkipRows() counts the number of header lines in pep.txt and psm.txt files
### returns value to be used for skiprows argument in pd.read_csv()
def getSkipRows(inFile):
    with open(inFile, "r") as f:
        lines = f.readlines()
        headerLines = 0
        for l in lines:
            if l.startswith("#"):
                headerLines += 1
    return headerLines

skipPepRows = getSkipRows(args.pepTxtFile)
skipPsmRows = getSkipRows(args.psmTxtFile)


### Read pep.txt and psm.txt files into pandas df & skip header lines
pepTxtDf = pd.read_csv(args.pepTxtFile, skiprows=skipPepRows, sep="\t")
psmTxtDf = pd.read_csv(args.psmTxtFile, skiprows=skipPsmRows, sep="\t")
lblPCTdf = pd.read_csv(args.lblPCT, sep="\t")

### Convert filtering threshold
pepAndPSMThresh = float(args.threshold) * 1000
proThresh = float(args.threshold)

### Parse psm.txt file and count total PSMs and number of labeled PSMs
countLabPSMs = 0
countTotalPSMs = 0
for psmName, psm in psmTxtDf.iloc[:, np.r_[-3, 8]].values:
    splitPsmName = psm.split("_")
    stripValue = splitPsmName[1].rstrip("Pct")
    if not psmName.startswith("{Rev_"):
        countTotalPSMs += 1
        if int(stripValue) >= pepAndPSMThresh:
            countLabPSMs += 1

### Parse pep.txt file and count total peptides and number of labeled peptides
countTotalPeptides = 0
countLabPeptides = 0
for pepName, pep in pepTxtDf.iloc[:, np.r_[3, -1]].values:
    splitPepName = pep.rstrip("}").lstrip("{").split(",")
    pepsSubList = []
    countTotalPeptides += 1
    for item in splitPepName:
        stripItem = item.lstrip("'C13").rstrip("Pct").lstrip("_")
        pepsSubList.append(int(stripItem))
    if not pepName.startswith("{Rev_"):
        if any(x >= pepAndPSMThresh for x in pepsSubList):
            countLabPeptides += 1 
            
## Parse pro.txt file and count total # proteins
with open(args.proTxtFile, 'r') as inHandle:
    for line in inHandle:
       if line.startswith('#	Total_Proteins_After_Filtering ='):
        splitLine = line.rstrip('\n').split(' ')
        countTotalPros = int(splitLine[-1])

## Parse LabelPCTCount file and count # labeled proteins
countLabPros = 0
for proName, enrichment in zip(lblPCTdf.iloc[:, 0].values, lblPCTdf.iloc[:, 1:-1].values.tolist()):
    splitName = proName.split(',')
    splitEnrich = [x.split(',') for x in enrichment if isinstance(x, str)]
    floatsList = [float(j) for k in splitEnrich for j in k]
    enrichArray = np.array(floatsList).flatten()
    if not np.all(enrichArray == 1.07):
        if any(y >= proThresh for y in enrichArray):
            countLabPros += 1


### Print summary stats
print("total PSMs = " + str(countTotalPSMs))
print("labeled PSMs = " + str(countLabPSMs))
print("total peptides = " + str(countTotalPeptides))
print("labeled peptides = " + str(countLabPeptides))
print("total proteins = " + str(countTotalPros))
print("labeled proteins = " + str(countLabPros))