#!/usr/bin/env python3
import pandas as pd
import numpy as np
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


def parseM8(inFile):
    m8df = pd.read_csv(inFile, sep = '\t', usecols = [0, 1, 13])

    taxLst = []
    for prot, cont, taxon in m8df.itertuples(index=False):
        if not pd.isna(taxon):  
            splitTax = taxon.split(';')
            taxLst.append(splitTax[0])

    fakeTaxList = list(range(len(taxLst)))
    fakeCountList = list(range(len(taxLst)))

    return taxLst, fakeTaxList, fakeCountList

taxList, fakeTaxonomy, fakeCount = parseM8(args.inFile)

with open(args.outFile, 'w') as outHandle:
    for t, fT, fC in zip(taxList, fakeTaxonomy, fakeCount):
        outHandle.write('{}\t{}\t{}\n'.format(fT, t, fC))
