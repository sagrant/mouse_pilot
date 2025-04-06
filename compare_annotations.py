#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

"""
Generate a data frame that compares all of the annotations generated with different softwares with one another, 
do they agree?

Inputs:
- nr: mmseqs output
- c: consensus_annotations_contigs.py output 
- g: consensus_annotations_genes.py output 
- b: dictionary that links bins to contigs
- s: gtdbtk output

Output:
- o: complete comparison data frame 
"""

class generateDictionaries():

    """Generate required lookup dictionaries

    The assembledge of the output dataframe relies on looking up all taxonomic annotations based on their 
    assoicated contig IDs with the dict.get() method.

    ...

    Attributes
    ----------
    mmseqsData : DataFrame
        data frame including contig ID, annotation rank, consensus annotation at the contig level, and complete ranks from mmseqs output
    binData : DataFrame
        data frame including bin name, complete ranks, and annotation at the bin level, among other fields from gtdbkt output 
    genesData : DataFrame
        data frame including gene ID, complete ranks, and consensus annotation at the gene level 
    consensusContigs : DataFrame
        data frame including contig ID, gene ID, consensus annotation at the contig level, annotation rank, and complete ranks
    """

    def __init__(self, mmseqsData, binData, genesData, consensusContigs):
        self.mmseqsData = mmseqsData
        self.binData = binData
        self.genesData = genesData
        self.consensusContigs = consensusContigs

    def parseMMseqsData(self):
        """ Parse mmseqs output. Returns two dictionaries
        
        Parameters
        ----------
        mmseqsData : DataFrame
            data frame including contig ID, annotation rank, consensus annotation at the contig level, and complete ranks from mmseqs output
    
        Returns 
        ----------
        getMMseqsRanksDict : dict
            Dictionary to be used to retrieve complete annotation ranks
        getMMseqsAnnotDict : dict
            Dictionary to be used to retrieve mmseqs contig level annotation 
        """
        getMMseqsRanksDict = {}
        getMMseqsAnnotDict = {}
        for mmContig, mmAnnot, mmRanks in self.mmseqsData.iloc[:, np.r_[0, 2, 3]].itertuples(index=False):
            mmseqs = []
            if isinstance(mmRanks, str):
                splitRnks = mmRanks.split(";")
            getMMseqsRanksDict[mmContig] = splitRnks
            getMMseqsAnnotDict[mmContig] = mmAnnot
        return getMMseqsRanksDict, getMMseqsAnnotDict
    
    def parseBinData(self):
        """ Parse gtdbtk output. Returns dictionary
        
        Parameters
        ----------
        binData : DataFrame
            data frame including bin name, complete ranks, and annotation at the bin level, among other fields from gtdbkt output 
    
        Returns 
        ----------
        binnedTaxaDict : dict
            Dictionary to be used to retrieve bin level annotation 
        """
        binnedTaxaDict = {}
        for genomeBin, classification in self.binData.iloc[:, 0:2].itertuples(index=False):
            binID = genomeBin.split(".")[-1]
            classificationList = classification.split(";")
            binnedTaxaDict[binID] = classificationList
        return binnedTaxaDict

    def parseGenesData(self):
        """ Parse consensus annotations genes output. Returns several dictionaries
        
        Parameters
        ----------
        genesData : DataFrame
            data frame including gene ID, complete ranks, and consensus annotation at the gene level 
    
        Returns 
        ----------
        genes2contigsDict : dict
            Dictionary to be used to retrieve all gene IDs associated with each contig
        geneRanksDict : dict
            Dictionary to be used to retrieve all gene level annotation complete based on gene ID
        geneAnnotsDict : dict
            Dictionary to be used to retrieve tuple with (annotation, taxonomic rank) for each gene 
        """
        self.genesData["Contig"] = self.genesData["Gene"].apply(lambda row: "_".join(row.rsplit("_")[0:3]))
        genes2contigsDict = defaultdict(list)
        geneRanksDict = {}
        geneAnnotsDict = {}
        for gene, geneRanks, annot, contigName in self.genesData.itertuples(index=False):
            genes2contigsDict[contigName].append(gene)
            geneRanksDict[gene] = geneRanks
            geneAnnotsDict[gene] = annot
        return genes2contigsDict, geneRanksDict, geneAnnotsDict

    def parseContigData(self, mmseqsRanksDict, mmseqsAnnotationsDict):
        allRanksDict = {}
        contigAnnotsDict = {}
        for conContig, conAnnot, rankLevel, conRanks in self.consensusContigs.iloc[:, np.r_[0, 2, 3, 4]].itertuples(index=False):
            customRanks = conRanks.split(",")
            custom = []
            for x in customRanks:
                trimx = x.rstrip(']"]').lstrip(" '['").rstrip("'")
                custom.append(trimx)
            mmseqsRnks = mmseqsRanksDict.get(conContig)
            mmseqsAnnot = mmseqsAnnotationsDict.get(conContig)
            allRanksDict[conContig] = [custom, mmseqsRnks]
            contigAnnotsDict[conContig] = [conAnnot, mmseqsAnnot]
        return allRanksDict, contigAnnotsDict

class makeComparison():

    def __init__(self, bins2contigsDictionary, binnedTaxaDict, ranksDict, contigAnnotsDict, genes2contigsDict, geneRanksDict, geneAnnotsDict):
        self.bins2contigsDictionary = bins2contigsDictionary
        self.binnedTaxaDict = binnedTaxaDict
        self.ranksDict = ranksDict
        self.contigAnnotsDict = contigAnnotsDict
        self.genes2contigsDict = genes2contigsDict
        self.geneRanksDict = geneRanksDict
        self.geneAnnotsDict = geneAnnotsDict

    def getGeneLevel(self, genes):
        for gene in genes:
            geneLevelRanks = self.geneRanksDict.get(gene)
            geneLevelAnnots = self.geneAnnotsDict.get(gene)
        return geneLevelRanks, geneLevelAnnots, gene

    def compareAll(self):
        annotsList = []
        for binName, binContigs in self.bins2contigsDictionary.items():
            binIdentifier = binName.replace(".", "_").split("_")[1]
            binLevelAnnots = self.binnedTaxaDict.get(binIdentifier)
            for contig in binContigs:
                contigLevelRanks = self.ranksDict.get(contig)
                contigLevelAnnots = self.contigAnnotsDict.get(contig)
                geneIDs = self.genes2contigsDict.get(contig)
                if geneIDs and contigLevelRanks:
                    geneRanks, geneTaxon, geneName = self.getGeneLevel(geneIDs)
                    if isinstance(geneTaxon, str):
                        geneAnnotation = geneTaxon.lstrip("('").split(",")[0].rstrip("'")
                        annotsList.append([binName, binLevelAnnots, contig, contigLevelRanks[1], contigLevelAnnots[1], contigLevelRanks[0], contigLevelAnnots[0], geneName, geneAnnotation, geneRanks])
        annotationsDf = pd.DataFrame(annotsList).rename(columns={0: "Bin",1: "Bin_Level_Ranks",2: "Contig",3: "mmeqs_Contig_Ranks",4: "mmseqs_Annotation",5: "Custom_Contig_Ranks",6: "Custom_Contig_Annotation",7: "Gene",8: "Gene_Level_Annotation",9: "Gene_Level_Ranks"})
        sortedAnnotsDf = annotationsDf.sort_values(by=["Contig"])
        return sortedAnnotsDf

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-nr", "--nrFile")
    parser.add_argument("-c", "--contigAnnotations")
    parser.add_argument("-g", "--geneAnnotations")
    parser.add_argument("-b", "--bcDict")
    parser.add_argument("-s", "--gtdbtk")
    parser.add_argument("-o", "--outFile")
    args = parser.parse_args()

    mmseqsDf = pd.read_csv(args.nrFile, sep="\t", header=None, usecols=[0, 2, 3, 8], names=["Contig", "Annotation_Rank", "Taxon", "Ranks"])

    consensusGenesDf = pd.read_csv(args.geneAnnotations, sep=",", header=0, names=["Gene", "Ranks", "Consensus"])

    consensusContigsDf = pd.read_csv(
        args.contigAnnotations, sep=",", header=0, names=["Contig", "Gene", "Consensus", "Consensus_Rank", "All_Ranks"])

    bins2contigs = pd.read_csv(args.bcDict, dtype=str).to_dict(orient="list")
    bins2contigsDict = {key: [x for x in value if not pd.isna(x)] for key, value in bins2contigs.items()}

    binsDf = pd.read_csv(args.gtdbtk, header=0, sep="\t")

    makeDicts = generateDictionaries(mmseqsDf, binsDf, consensusGenesDf, consensusContigsDf)
    mmRanks, mmAnnotations = makeDicts.parseMMseqsData()
    binTaxa = makeDicts.parseBinData()
    genes2contigs, geneRanks, geneAnnots = makeDicts.parseGenesData()
    allRanks, contigAnnots = makeDicts.parseContigData(mmRanks, mmAnnotations)

    comparison = makeComparison(bins2contigsDict, binTaxa, allRanks, contigAnnots, genes2contigs, geneRanks, geneAnnots)
    completeComparisonDf = comparison.compareAll()
    completeComparisonDf.to_csv(args.outFile, index = False)

if __name__ == "__main__":
    main()

