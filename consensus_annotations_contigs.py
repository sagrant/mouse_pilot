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
Assign consensus annotations to contigs and proteins based on gb_taxonomy_tools and DIAMOND outputs 

args
-i = gene annotations file for each sample
-o = consensus annotations & ranks to be input for compare_annotations.py
-m = missed annotations
"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inFile")
parser.add_argument("-o", "--outFile")
parser.add_argument("-m", "--missed")
args = parser.parse_args()


class findConsensusAnnotations_contigs():

    def __init__(self, inDataFrame):
        self.inDataFrame = inDataFrame

    
    def checkConsensusMatch(self):
        ### Get contig IDs from gene IDs and groupby contig
        ### Create dictionary with contig ID as keys and consensus annotations as values
        self.inDataFrame["Contig"] = self.inDataFrame.iloc[:, 0].apply(lambda row: "_".join(row.rsplit("_")[0:3]))
        gbContig = self.inDataFrame.groupby(["Contig"])
        
        ranksDict = {}  # use to get ranks if required to step back
        taxonHashDict = {}  # will become lookUpDict
        mainHashList = []  # will become hashLists
        noConsensusContigs = set()
        geneConsensus = []
        writeRanksOutDict = {}
        writeGenesOutDict = {}
        for x, group in gbContig:
            geneLevelConsensus = []
            rankHashes = []
            writeGenesOutDict[x] = group["Gene"].values.tolist()
            for gene, rnks, consensus, contig in group.itertuples(index=False):
                splitRanks = rnks.split("|")
                nHash = hash("nan")

                if not isinstance(consensus, float):
                    if len(splitRanks) > 1:
                        geneLevelConsensus.append(consensus)
                        hashedRanks = []
                        indexedRankList = [0, 0, 0, 0, 0, 0, 0]
                        indexedRankList[0] = splitRanks[0]
                        indexedRankList[1] = splitRanks[4]
                        indexedRankList[2] = splitRanks[7]
                        indexedRankList[3] = splitRanks[10]
                        indexedRankList[4] = splitRanks[13]
                        indexedRankList[5] = splitRanks[17]
                        indexedRankList[6] = splitRanks[19]
                        indexedRankList.reverse()

                        for taxon in indexedRankList:
                            hashedTaxon = hash(taxon)
                            hashedRanks.append(hashedTaxon)
                            taxonHashDict[x, hashedTaxon] = taxon
                            writeRanksOutDict[x, hashedTaxon] = [indexedRankList, taxon]
                        rankHashes.append(hashedRanks)
                    else:
                        rankHashes.append(nHash)
                else:
                    noConsensusContigs.add(x)
                ranksDict[x] = rankHashes
            mainHashList.append([x, rankHashes, geneLevelConsensus])
            geneConsensus.append([x, geneLevelConsensus])
        return ranksDict, nHash, taxonHashDict, mainHashList, noConsensusContigs, geneConsensus, writeRanksOutDict, writeGenesOutDict

    ### Account for genera/families/orders/classes/phyla where # NaNs surpasses threshold
    def eliminateNAs(self, someRank, contig, nanHashed):
        ### the 0th index of each item in misses is 0, which is used to disqualify any set ranks where NaN was majority
        misses = [0]
        if nanHashed in someRank:
            countNaNs = someRank.count(nanHashed)
            if (countNaNs / len(someRank)) > 0.5:
                misses.append(contig)
                return misses
            else:
                return someRank
        else:
            return someRank


    ### Function returns majority taxon for genera, family, order...
    def findMajorityConsensus(self, ranks, contig):
        if ranks and ranks[0] != 0:
            x = Counter(ranks)
            majority = x.most_common()
            return majority[0]
        else:
            return 0, 0


    ### Function returns majority taxon for genera, family, order...
    def findMajority(self, ranks, contig):
        ### when ranks[0] = 0 that signifies that majority of that rank was NaN
        if ranks and ranks[0] != 0:
            x = Counter(ranks)
            majority = x.most_common()
            return majority[0]
        else:
            return 0, 0


    ### Evaluate findMajoirty: use threshold to determine which taxon is majority
    def evaluateMajorityConsensus(self, consTaxon, consensusTaxonList, cont, countMaj, numAnnotations):
        confValDict = {}
        consensusAgree = False
        if consTaxon != 0:
            splitConsTax = consTaxon.split("', ")[0].lstrip("('")
            splitConsLevel = consTaxon.split(", ")[1].rstrip("')")  ## BU3 CC: IndexError: list index out of range
            if (countMaj / numAnnotations) > 0.5:
                consensusTaxonList.append([(cont, hash(splitConsTax)), splitConsLevel])
                consensusAgree = True
                confValDict[cont] = countMaj / numAnnotations
        return consensusAgree


    ### Evaluate findMajoirty: use threshold to determine which taxon is majority
    def evaluateMajority(self, taxonList, cont, emptyList, numAnnotations, keyWord):
        confValDict = {}
        consensusFound = False
        currentMajority, currentCount = self.findMajority(taxonList, cont)
        ### when current majority = 0, the majority hash was NaN
        if currentMajority != 0 and numAnnotations != 0:
            if (currentCount / numAnnotations) > 0.5:
                emptyList.append([(cont, currentMajority), keyWord])
                consensusFound = True
                confValDict[cont] = currentCount / numAnnotations
        return consensusFound, confValDict

    def findConsensus(self, hashLists, ranksDictionary, NaNHash):
        ### Loop through hashed taxonomic ranks and determine majority
        consensusSpecies = []
        consensusGenera = []
        consensusFamily = []
        consensusOrder = []
        consensusClass = []
        consensusPhylum = []
        consensusKingdom = []
        missedList = []
        checkConsensusAnnots = []
        confidenceDict = {}
        for c, hList, consensusList in hashLists:
            species = []
            genera = []
            families = []
            orders = []
            clsses = []
            phyla = []
            domains = []
            getRanks = []
            totalLen = len(hList)
            consenusTaxon, countTaxon = self.findMajorityConsensus(consensusList, c)
            consensusBool = self.evaluateMajorityConsensus(consenusTaxon, checkConsensusAnnots, c, countTaxon, totalLen)
            if not consensusBool:
                for contigGroup in hList:
                    ranks = ranksDictionary.get(c)
                    getRanks.append(ranks)
                for subList in getRanks:
                    for s in subList:
                        if s != NaNHash and len(s) > 0:
                            species.append(s[0])
                        if s != NaNHash and len(s) > 1:
                            genera.append(s[1])
                        if s != NaNHash and len(s) > 2:
                            families.append(s[2])
                        if s != NaNHash and len(s) > 3:
                            orders.append(s[3])
                        if s != NaNHash and len(s) > 4:
                            clsses.append(s[4])
                        if s != NaNHash and len(s) > 5:
                            phyla.append(s[5])
                        if s != NaNHash and len(s) > 6:
                            domains.append(s[6])

                filteredSpecies = self.eliminateNAs(species, c, NaNHash)
                filteredGenera = self.eliminateNAs(genera, c, NaNHash)
                filteredFamilies = self.eliminateNAs(families, c, NaNHash)
                filteredOrders = self.eliminateNAs(orders, c, NaNHash)
                filteredClasses = self.eliminateNAs(clsses, c, NaNHash)
                filteredPhyla = self.eliminateNAs(phyla, c, NaNHash)
                filteredDomains = self.eliminateNAs(domains, c, NaNHash)

                speciesBool, speciesConfVals = self.evaluateMajority(filteredSpecies, c, consensusSpecies, totalLen, "Species")
                confidenceDict.update(speciesConfVals)
                if speciesBool:
                    continue

                generaBool, generaConfVals = self.evaluateMajority(filteredGenera, c, consensusGenera, totalLen, "Genus")
                confidenceDict.update(generaConfVals)
                if generaBool:
                    continue
                familyBool, familyConfVals = self.evaluateMajority(filteredFamilies, c, consensusFamily, totalLen, "Family")
                confidenceDict.update(familyConfVals)
                if familyBool:
                    continue
                orderBool, orderConfVals = self.evaluateMajority(filteredOrders, c, consensusOrder, totalLen, "Order")
                confidenceDict.update(orderConfVals)
                if orderBool:
                    continue
                classBool, classConfVals = self.evaluateMajority(filteredClasses, c, consensusClass, totalLen, "Class")
                confidenceDict.update(classConfVals)
                if classBool:
                    continue
                phylaBool, phylaConfVals = self.evaluateMajority(filteredPhyla, c, consensusPhylum, totalLen, "Phylum")
                confidenceDict.update(phylaConfVals)
                if phylaBool:
                    continue
                kingdomsBool, kingdomsConfVals = self.evaluateMajority(filteredDomains, c, consensusKingdom, totalLen, "Kingdom")
                confidenceDict.update(kingdomsConfVals)
                if kingdomsBool:
                    continue
                else:
                    missedList.append(c)

        concatList = (
            consensusSpecies
            + consensusGenera
            + consensusFamily
            + consensusOrder
            + consensusClass
            + consensusPhylum
            + consensusKingdom
            + checkConsensusAnnots
            )
        return concatList, missedList, confidenceDict


    def unHash(self, allConcatList, lookupDict):
        ### Get taxon name and assoicated contig/protein based on hash value
        findDict = {}
        for item in allConcatList:
            res = lookupDict.get(item[0])
            findDict[item[0][0]] = [res, item[1]]
        return findDict


    def calculateStats(self, findDict):
        speciesCount = 0
        genusCount = 0
        familyCount = 0
        orderCount = 0
        classCount = 0
        phylumCount = 0
        kingdomCount = 0
        for values in findDict.values():
            stripVals = values[1].lstrip("'")
            if stripVals == "Species":
                speciesCount += 1
            if stripVals == "Genus":
                genusCount += 1
            if stripVals == "Family":
                familyCount += 1
            if stripVals == "Order":
                orderCount += 1
            if stripVals == "Class":
                classCount += 1
            if stripVals == "Phylum":
                phylumCount += 1
            if stripVals == "Kingdom":
                kingdomCount += 1
        
        percentS = spCnt / len(findDict) * 100
        percentG = genCnt / len(findDict) * 100
        percentF = famCnt / len(findDict) * 100
        percentO = ordCnt / len(findDict) * 100
        percentC = clsCnt / len(findDict) * 100
        percentP = phyCnt / len(findDict) * 100
        percentK = kingCnt / len(findDict) * 100
        return percentS, percentG, percentF, percentO, percentC, percentP, percentK

def main():
    inDataFrame = pd.read_csv(
    args.inFile, sep=",", header=0, dtype=str, names=["Gene", "Ranks", "Consensus"]
)   
    finder_contigs = findConsensusAnnotations_contigs(inDataFrame)
    ranksDict, nanHash, LUDict, hashLists, noCons, geneConsensusList, ranksOutDictionary, genesOut = finder_contigs.checkConsensusMatch()
    allAnnotsList, missedAnnotsList, confidenceStats = finder_contigs.findConsensus(hashLists, ranksDict, nanHash)
    foundTaxaDict = finder_contigs.unHash(allAnnotsList, LUDict)
    pctS, pctG, pctF, pctO, pctC, pctP, pctK = finder_contigs.calculateStats(foundTaxaDict)

    ### Print summary stats to STDOUT
    print(args.inFile)
    print("total # contigs in = " + str(len(ranksOutDictionary)))
    print("total # assigned annotations = " + str(len(foundTaxaDict)))
    print("total # missed annotations = " + str(len(missedAnnotsList)))
    print("percent assigned = " + str((len(foundTaxaDict) / len(ranksOutDictionary) * 100)))
    print(" ")

    print("percent species = " + str(pctS))
    print("percent genus = " + str(pctG))
    print("percent family = " + str(pctF))
    print("percent order = " + str(pctO))
    print("percent class = " + str(pctC))
    print("percent phylum = " + str(pctP))
    print("percent kingdom = " + str(pctK))
    print(" ")

    with open(args.missed, "w") as missedHandle:
        for item in missedAnnotsList:
            missedHandle.write(str(item) + "\n")

    outList = []
    for contigEntry in allAnnotsList:
        getRanksOut = ranksOutDictionary.get(contigEntry[0])
        getTaxonOut = LUDict.get(contigEntry[0])
        getGeneOut = genesOut.get(contigEntry[0][0])
        outList.append([contigEntry[0][0], getGeneOut[0], getTaxonOut, contigEntry[1].lstrip("'"), getRanksOut[0]])

    outDf = pd.DataFrame(outList).rename(columns={0: "Contig", 1: "Genes", 2: "Consensus_Annotation", 3: "Level", 4: "Ranks"})
    explodeOut = outDf.explode("Genes")
    explodeOut.to_csv(args.outFile, index=False)

if __name__ == "__main__":
    main()