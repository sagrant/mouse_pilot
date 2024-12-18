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
-o = output file to be used as input for consensus_annotations_genes.py (.txt)

'''

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gbFile')
parser.add_argument('-m', '--m8File')
parser.add_argument('-o', '--outFile')
args = parser.parse_args()

class generateMerged():

    def __init__(self, m8_df, gb_df):
        self.m8Df = m8_df
        self.gbDf = gb_df


    def makeRanksDict(self):
        ### iterate through gb_taxonomy_tools output and create dictionary where keys are taxon ID and values are all taxonomic ranks (including NaNs)
        taxDict = {}
        for val1, taxid, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16, val17, val18, val19, val20, val21, val22, val23, val24, val25 in self.gbDf.itertuples(index=False):
            taxDict[taxid] = str(val4) + '|' + \
                str(val5) + '|' + str(val6) + '|' + str(val7) + '|' + str(val8) + '|' + str(val9) + '|' + str(val10) + '|' + \
                str(val11) + '|' + str(val12) + '|' + str(val13) + '|' + str(val14) + '|' + str(val15) + '|' + str(val16) + '|' + str(val17) + '|' + \
                str(val18) + '|' + str(val19) + '|' + str(val20) + '|' + str(val21) + '|' + str(val22) + '|' + str(val23) + '|' + str(val24) + '|' + str(val25)
        return taxDict


    def getIDs(self, taxonDict):
        ### get contig and protein ID from m8 based on taxon IDs
        lookupList = []
        countMatches = 0
        countNoMatches = 0
        countNoID = 0
        for cont, prot, tax in self.m8Df.itertuples(index = False):
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


def main():
    m8df = pd.read_csv(args.m8File, sep='\t', header=None, usecols=[0, 1, 13], dtype=str).rename(columns={0: 'Contig', 1: 'Protein', 2: 'Taxon'})
    gbDf = pd.read_csv(args.gbFile, sep = '\t', header = None, na_values = 'n', dtype = str).rename(columns = {0: 'Taxon'})

    merger = generateMerged(m8df, gbDf)
    taxonomyDict = merger.makeRanksDict()
    luList, cntMatch, cntNoMatch, countNoid = merger.getIDs(taxonomyDict)

    ### Print summary stats to STDOUT
    print(f'Sample = ' + str(args.m8File.split('_')[0])) ### Assumes sample name is same as file name prefix
    print('number matched = ' + str(cntMatch))
    print('number not matched = ' + str(cntNoMatch))
    print('number no ID in m8 = ' + str(countNoid))

    ### Write output file
    m8df['Ranks'] = luList
    m8df.to_csv(args.outFile, index = False)

if __name__ == "__main__":
    main()