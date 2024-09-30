#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

'''
02-26-2024
Merge m8 output with gb_taxonomy_tools output on Taxon ID
Create output file to be input for consensus_annotations_genes.py

args:
-g = gb_taxonomy_tools output file (.taxonomy)
-m = DIAMOND output (.m8)
-o = output file to be used as input for consensus_annotations_genes.py

'''

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gbFile')
parser.add_argument('-m', '--m8file')
parser.add_argument('-o', '--outFile')
args = parser.parse_args()

def readInFiles(m8File, gbFile):
    m8df = pd.read_csv(m8File, sep='\t', header=None, usecols=[0, 1, 13], dtype=str).rename(columns={0: 'Contig', 1: 'Protein', 2: 'Taxon'})
    gbDf = pd.read_csv(gbFile, sep = '\t', header = None, na_values = 'n', dtype = str).rename(columns = {0: 'Taxon'})
    return m8df, gbDf

m8_df, gb_df = readInFiles(args.m8file, args.gbFile)

def makeRanksDict(gbData):
    ### iterate through gb_taxonomy_tools output and create dictionary where keys are taxon ID and values are all taxonomic ranks (including NaNs)
    taxDict = {}
    for val1, taxid, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16, val17, val18, val19, val20, val21, val22, val23, val24, val25 in gbData.itertuples(index=False):
        taxDict[taxid] = str(val4) + '|' + \
            str(val5) + '|' + str(val6) + '|' + str(val7) + '|' + str(val8) + '|' + str(val9) + '|' + str(val10) + '|' + \
            str(val11) + '|' + str(val12) + '|' + str(val13) + '|' + str(val14) + '|' + str(val15) + '|' + str(val16) + '|' + str(val17) + '|' + \
            str(val18) + '|' + str(val19) + '|' + str(val20) + '|' + str(val21) + '|' + str(val22) + '|' + str(val23) + '|' + str(val24) + '|' + str(val25)
    return taxDict

taxonDict = makeRanksDict(gb_df)

def getIDs(m8df):
    ### get contig and protein ID from m8 based on taxon IDs
    lookupList = []
    countMatches = 0
    countNoMatches = 0
    countNoID = 0
    for cont, prot, tax in m8df.itertuples(index = False):
        if not pd.isna(tax):
            splitTax = tax.split(';')
            res = taxonDict.get(splitTax[0])
            if res != None:
                countMatches += 1
                lookupList.append(res)
            if res == None:
                countNoMatches += 1
                lookupList.append('No Match') ### signifies script couldnt find a match in taxon IDs between gb_taxonomy_tools output and m8 output
        else:
            lookupList.append('No ID in m8') ### signifies no taxon ID was present in m8
            countNoID += 1
    return lookupList, countMatches, countNoMatches, countNoID

lookUpData, countMatch, countNoMatch, countNoID = getIDs(m8_df)

### Print summary stats to STDOUT
print(f'Sample = ' + str(args.m8file.split('_')[0])) ### Assumes sample name is same as file name
print('number matched = ' + str(countMatch))
print('number not matched = ' + str(countNoMatch))
print('number no ID in m8 = ' + str(countNoID))

### Write output file
m8_df['Ranks'] = lookUpData
m8_df.to_csv(args.outFile, index = False)
