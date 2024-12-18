#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from collections import Counter
import os
import sys

### Set hash seed so results are reproducible
hashseed = os.getenv("PYTHONHASHSEED")
if not hashseed:
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

"""
02-26-2024
Assign consensus annotations to contigs and proteins based on gb_taxonomy_tools and DIAMOND outputs 

args
-i = input file is output from merge_taxon_IDs.py (merged file)
-r = output csv file with all assigned annotations and taxonomic ranks
-m = output csv with list of genes that received no consensus annotation

"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inFile")
parser.add_argument("-r", "--ranks")
parser.add_argument("-m", "--missed")
args = parser.parse_args()


class findConsensusAnnotations_genes():

    def __init__(self, infile):
        self.infile = infile

    def parseInData(self):
        ### Parse input file
        ### Return groupby object and main dataframe with all data
        inDataFrame = pd.read_csv(self.infile, sep=",", header=0, dtype=str, names=["Gene", "Protein", "Taxon", "Ranks"]).fillna(0)

        ### Groupby gene to get gene-level annotations
        gbGene = inDataFrame.groupby("Gene")
        return gbGene, inDataFrame


    def getTaxonomyDict(self, df, colName):
        ### Iterate through groupby object to save taxonomic annotation predictions 
        ### Return dictionary with gene ID as keys and taxonomic ranks as values
        ### Return list of genes that had no prediction 
        annotDict = {}
        emptyGenes = []
        for x, group in df:
            ### if all annotations for a contig are all NaN, don't try to find consensus
            if not (group["Taxon"] == 0).all():
                gene = group[colName]
                taxa = group["Ranks"].str.split("|").values.tolist()
                annotDict[gene.iloc[0]] = taxa
            else:
                emptyGenes.append(group[colName].values[0])
        return annotDict, emptyGenes


    def hashTaxa(self, geneTaxonomy):
        ### Hash taxa in gene dictionary
        ### Return hashed taxonomic names
        ### Return lookup dictionary so taxon names can be retreived based on hash value
        ### Return hash value for 'nan' to be used in eliminateNAs()
        hashLists = []
        lookUpDict = {}
        for k, v in geneTaxonomy.items():
            tempHashList = []
            for rankList in v:
                nanHash = hash("nan")
                if len(rankList) > 1:
                    hashedRankList = []
                    ### indexed taxonlist contains species genus, family, order, class, and phylum 
                    indexedTaxonList = [0, 0, 0, 0, 0, 0, 0]
                    indexedTaxonList[0] = rankList[0]
                    indexedTaxonList[1] = rankList[4]
                    indexedTaxonList[2] = rankList[7]
                    indexedTaxonList[3] = rankList[10]
                    indexedTaxonList[4] = rankList[13]
                    indexedTaxonList[5] = rankList[17]
                    indexedTaxonList[6] = rankList[19]
                    indexedTaxonList.reverse()
                    for taxon in indexedTaxonList:
                        hashItems = hash(taxon)
                        hashedRankList.append(hashItems)
                        countNaNs = hashedRankList.count(nanHash)
                        lookUpDict[k, hashItems] = taxon
                    tempHashList.append(hashedRankList)
                else:
                    tempHashList.append([nanHash])
            hashLists.append([k, tempHashList])
        return hashLists, lookUpDict, countNaNs


    def eliminateNAs(self, someRank, g, cnan):
        ### Account for genera/families/orders/classes/phyla where # NaNs surpasses threshold
        ### Return ranks if # NaNs do NOT surpass threshold, return 0 when # NaNs DO surpass threshold
        misses = [0]
        if cnan in someRank:
            countNaNs = someRank.count(cnan)
            if (countNaNs / len(someRank)) > 0.5:
                misses.append(g)
                return misses
            else:
                return someRank
        else:
            ### the 0th index of each item in misses is 0, which is used to disqualify any set ranks where NaN was majority
            return someRank


    def findMajority(self, ranks, g):
        ### Find consensus annotation based on majority rule
        ### Function returns majority taxon for genera, family, order...
        ### when ranks[0] = 0 that signifies that majority of that rank was NaN
        if ranks[0] != 0:
            x = Counter(ranks)
            majority = x.most_common()
            return majority[0]
        else:
            return 0, 0


    def evaluateMajority(self, taxonList, consensusTaxonList, gene, numAnnotations, keyWord):
        ### Evaluate findMajoirty: use threshold to determine which taxon is majority
        confValDict = {}
        consensusFound = False
        currentMajority, currentCount = self.findMajority(taxonList, gene)
        ### when current majority = 0, the majority hash was NaN
        if currentMajority != 0:
            if (currentCount / numAnnotations) > 0.5:
                consensusTaxonList.append([(gene, currentMajority), keyWord])
                consensusFound = True
                confValDict[gene] = currentCount / numAnnotations
        return consensusFound, confValDict


    def findConsensus(self, hashList, cntNaNs):
        ### Loop through hashed taxonomic ranks and determine majority
        consensusSpecies = []
        consensusGenera = []
        consensusFamily = []
        consensusOrder = []
        consensusClass = []
        consensusPhylum = []
        consensusKingdom = []
        missedList = []
        confidenceDict = {}
        for gene, hList in hashList:
            ### hashList contains 5 sublists with hashed taxonomic ranks per gene
            ### Set k=5 in DIAMOND command --> get 5 taxonomic predictions per gene 
            species = []
            genera = []
            families = []
            orders = []
            clsses = []
            phyla = []
            domains = []
            for geneGroup in hList:
                ### Sort taxon hash values by taxonomic rank
                if len(geneGroup) > 0:
                    species.append(geneGroup[0])
                if len(geneGroup) > 1:
                    genera.append(geneGroup[1])
                if len(geneGroup) > 2:
                    families.append(geneGroup[2])
                if len(geneGroup) > 3:
                    orders.append(geneGroup[3])
                if len(geneGroup) > 4:
                    clsses.append(geneGroup[4])
                if len(geneGroup) > 5:
                    phyla.append(geneGroup[5])
                if len(geneGroup) > 6:
                    domains.append(geneGroup[6])
            filteredSpecies = self.eliminateNAs(species, gene, cntNaNs)
            filteredGenera = self.eliminateNAs(genera, gene, cntNaNs)
            filteredFamilies = self.eliminateNAs(families, gene, cntNaNs)
            filteredOrders = self.eliminateNAs(orders, gene, cntNaNs)
            filteredClasses = self.eliminateNAs(clsses, gene, cntNaNs)
            filteredPhyla = self.eliminateNAs(phyla, gene, cntNaNs)
            filteredDomains = self.eliminateNAs(domains, gene, cntNaNs)
            totalNumAnnotations = len(hList)

            ### if species level consensus is found, save it. Otherwise, check if there is a genus level consensus, and so on... (LCA method)
            speciesBool, speciesConfVals = self.evaluateMajority(filteredSpecies, consensusSpecies, gene, totalNumAnnotations, "Species")
            confidenceDict.update(speciesConfVals)
            if speciesBool:
                continue
            generaBool, generaConfVals = self.evaluateMajority(filteredGenera, consensusGenera, gene, totalNumAnnotations, "Genus")
            confidenceDict.update(generaConfVals)
            if generaBool:
                continue
            familyBool, familyConfVals = self.evaluateMajority(filteredFamilies, consensusFamily, gene, totalNumAnnotations, "Family")
            confidenceDict.update(familyConfVals)
            if familyBool:
                continue
            orderBool, orderConfVals = self.evaluateMajority(filteredOrders, consensusOrder, gene, totalNumAnnotations, "Order")
            confidenceDict.update(orderConfVals)
            if orderBool:
                continue
            classBool, classConfVals = self.evaluateMajority(filteredClasses, consensusClass, gene, totalNumAnnotations, "Class")
            confidenceDict.update(orderConfVals)
            if classBool:
                continue
            phylaBool, phylaConfVals = self.evaluateMajority(filteredPhyla, consensusPhylum, gene, totalNumAnnotations, "Phylum")
            confidenceDict.update(phylaConfVals)
            if phylaBool:
                continue
            kingdomsBool, kingdomsConfVals = self.evaluateMajority(filteredDomains, consensusKingdom, gene, totalNumAnnotations, "Kingdom")
            confidenceDict.update(kingdomsConfVals)
            if kingdomsBool:
                continue
            else:
                missedList.append(gene)
        
        ### Concatenate results into one list
        concatList = (
            consensusSpecies + consensusGenera + consensusFamily + consensusOrder + consensusClass + consensusPhylum + consensusKingdom
        )

        ### Calculate number of annotations at each taxonomic rank to compute summary stats
        consSpLen = len(consensusSpecies)
        consGeLen = len(consensusGenera)
        consFaLen = len(consensusFamily)
        consOrLen = len(consensusOrder)
        consClLen = len(consensusClass)
        consPhLen = len(consensusPhylum)
        consKiLen = len(consensusKingdom)
        return concatList, missedList, consSpLen, consGeLen, consFaLen, consOrLen, consClLen, consPhLen, consKiLen

    def unHash(self, allAnnots, lookUp):
        ### Determine taxonomic consensus annotation for each gene based on hash value
        ### return dictionary where keys are gene ID and values are (taxon name, taxonomic rank) 
        findDict = {}
        for item in allAnnots:
            ### Link consensus annotation taxon name to taxon hash value 
            res = lookUp.get(item[0])
            findDict[item[0][0]] = (res, item[1])
        findDf = pd.DataFrame().from_dict(findDict, orient="index").reset_index().rename(columns={"index": "Gene", 0: "Annotation"})
        return findDf, findDict


    def writeOut(self, missedTaxaList, mainData, fDict, missedOutFileName, ranksOutFileName):
        ### Write file containing genes that got no consensus
        for item in missedTaxaList:
            row1 = mainData[mainData["Gene"] == item]
            gene1 = row1.iloc[:, 0]
            prot1 = row1.iloc[:, 1]
            tax1 = row1.iloc[:, 2]
            rnk1 = row1.iloc[:, 3]
            noConsensusGenes = pd.DataFrame({"Gene": gene1.values, "Taxon_IDs": tax1.values, "Ranks": rnk1.values})
            noConsensusGenes.to_csv(missedOutFileName)

        ### Get assigned gene IDs so annotations can be assigned to contigs. 
        ### Get ranks so that if there is no contig-level consensus at the species level, a consensus can be found among genera or families, etc.
        ranksOut = []
        for gene, ranks in mainData.iloc[:, np.r_[0, 3]].itertuples(index=False):
            result = fDict.get(gene)
            ranksOut.append({"Gene": gene, "Ranks": ranks, "Consensus": result})

        ranksOutDf = pd.DataFrame(ranksOut)
        ranksOutDf.to_csv(ranksOutFileName, index=False)



def main():
    finder_genes = findConsensusAnnotations_genes(args.inFile)
    groupGenesDf, mainDf = finder_genes.parseInData()
    geneTaxonomyDict, emptyGenesList = finder_genes.getTaxonomyDict(groupGenesDf, 'Gene')
    mainHashList, lookupDictionary, countNaN = finder_genes.hashTaxa(geneTaxonomyDict)
    allAnnotsList, missedTaxa, SpeciesLen, GenusLen, FamilyLen, OrderLen, ClassLen, PhylumLen, KingdomLen = finder_genes.findConsensus(mainHashList, countNaN)
    foundDf, foundDict = finder_genes.unHash(allAnnotsList, lookupDictionary)
    finder_genes.writeOut(missedTaxa, mainDf, foundDict, args.missed, args.ranks)


    ## Write summary stats to STDOUT
    print(f'{args.inFile} summary stats')
    print(f"Number of unique genes = {len(geneTaxonomyDict)}")
    print(f"Percent assigned = {(len(foundDict)/len(geneTaxonomyDict))}")
    print(f"Number of missed genes = {len(missedTaxa)}")
    print(" ")
    print(f"Percent species level annotations = {(SpeciesLen/len(foundDict))*100}")
    print(f"Percent genus level annotations = {(GenusLen/len(foundDict))*100}")
    print(f"Percent family level annotations = {(FamilyLen/len(foundDict))*100}")
    print(f"Percent order level annotations = {(OrderLen/len(foundDict))*100}")
    print(f"Percent class level annotations = {(ClassLen/len(foundDict))*100}")
    print(f"Percent phylum level annotations = {(PhylumLen/len(foundDict))*100}")
    print(f"Percent kingdom level annotations = {(KingdomLen/len(foundDict))*100}")
    print(" ")

if __name__ == "__main__":
    main()