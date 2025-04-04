# mouse_pilot
This repository contains code and instructions to carry out a bioinformatics workflow for gut microbiome analysis following a Proteomic Stable Isotope Probing experiment. This workflow was designed to process results from an experiment that investigated the effect of prebiotics on the mouse gut microbiome. 

Prebiotics are carbohydrates that are degraded by beneficial bacteria in the gut microbiome, and they improve host health by providing a nutrient source for these beneficial populations. We do not have a way to identify the specific microbial species that degrade prebiotics, which restricts our ability to understand how we could support those populations. To address this knowledge gap, we use a method called Proteomic Stable Isotope probing, where a prebiotic composed of stable isotope-labeled carbons is administered to live mice. The stable isotope will be assimilated by taxa that are actively degrading the prebiotic and will be subsequently incorporated into microbial proteins within the gut microbiome. Because the stable isotope of carbon is heavier than the naturally occuring isotope, proteins synthesized by prebiotic-degrading species are identifiable with mass spectrometry. The heavy protein sequences identified with mass spec are then matched to taxonomic sequnces from a paired metagenome, and a link can be made between the prebiotic substrate and the taxa responsible for degrading it. 

There are two sides of this analysis: 
1. Generation of metagenome annotations
3. Analysis of proteomics results

Overview of entire workflow:
![github_figure1](https://github.com/user-attachments/assets/c4fffb37-0260-4cef-92d7-ae8ddbc915d4)


A diagram of the metagenome annotations workflow is displayed below:
![github_figures](https://github.com/user-attachments/assets/e97ef049-2e15-441c-8fb7-6c8cb72d2be2)

The goal of this analysis is to assign taxonomic annotations to genes, contigs, and bins. Multiple softwares are used to increase our confidence that each annotation is correct. Our approach combines the Least Common Ancestor (LCA) and majority rule methods to assign annotations to genes, contigs, or bins.

Description of each step in the metagenome annotations workflow:

1. Generate annotations at the contig and gene level:
   - 1.a) Generate contig-level annotations with `mmseqs` software.
   - 1.b) Generate gene-level annotations with `DIAMOND` software.
2. Generate a dummy `.taxid` file for `gb_taxonomy_tools`.
     - Normal `.taxid` files have 3 columns: GenBank ID, taxon ID, and count.
     - A dummy `.taxid` simulates the first and third columns with irrelevant integers. The second column contains valid taxon IDs.
     - This allows the researcher to retrieve all taxonomic ranks for each `DIAMOND` annotation without carrying out the entire `gb_taxonomy_tools` workflow.
3. Run `gb_taxonomy_tools` with the dummy `.taxid` file to get taxonomic ranks for each `DIAMOND` annotation.
4. Merge `DIAMOND` and `gb_taxonomy_tools` output.
     - The output file contains 4 columns: Gene ID, Protein ID, NCBI Taxon ID, and taxonomic ranks.
5. Generate consensus annotations for each gene:
   - Two output files are created:
     - `-o` file: Gene ID, taxonomic ranks, and gene-level consensus annotation.
     - `-m` file: Gene ID, list of taxonomic IDs output by `DIAMOND` for that gene, and ranks output by `gb_taxonomy_tools` for each taxonomic ID.
   - Summary statistics printed to STDOUT include:
     - Number of unique genes.
     - Percent of genes that were assigned a consensus annotation.
     - Number of genes that were missed.
     - Breakdown of annotation levels (species, genus, family, etc.).
6. Generate consensus annotations for each contig:
   - Two output files are created:
     - `-o` file: Contig ID, Gene ID, contig-level consensus annotation, and taxonomic ranks.
     - `-m` file: Contig IDs not assigned a consensus annotation.
   - Summary statistics printed to STDOUT include:
     - Number of input contigs.
     - Number of contigs assigned a consensus annotation.
     - Number of contigs missed.
     - Breakdown of annotation levels (species, genus, family, etc.).
7. Combine all annotation data to facilitate comparison:
   - Determine if annotations output by different softwares agree with one another.
8. Generate consensus annotations for each bin:
   - If there is no consensus within a bin, consensus annotations are provided at the gene level.
