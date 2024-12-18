#!/usr/bin/env python3
import pandas as pd
import argparse

'''
02-14-2024
Generate dummy file for gb_taxonomy_tools

args:
-i = DIAMOND output (.m8 file)
-o = dummy file to use as input for gb_taxonomy_tools 
'''

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile')
parser.add_argument('-o', '--outFile')
args = parser.parse_args()

class getDummy():

    def __init__(self, m8Df):
        self.m8Df = m8Df


    def parseM8(self,):
        taxLst = []
        for prot, cont, taxon in self.m8Df.itertuples(index=False):
            if not pd.isna(taxon):  
                splitTax = taxon.split(';')
                taxLst.append(splitTax[0])

        fakeTaxList = list(range(len(taxLst)))
        fakeCountList = list(range(len(taxLst)))
        return taxLst, fakeTaxList, fakeCountList


def main():
    ### Read DIAMOND output 
    m8df = pd.read_csv(args.inFile, sep = '\t', usecols = [0, 1, 13])

    parseFunc = getDummy(m8df)
    taxonList, fakeTaxonList, fakeCountList = parseFunc.parseM8()

    #write concatenated paired end FASTQs to new file
    with open(args.outFile, 'w') as outHandle:
        for t, fT, fC in zip(taxonList, fakeTaxonList, fakeCountList):
            outHandle.write('{}\t{}\t{}\n'.format(fT, t, fC))

if __name__ == "__main__":
    main()