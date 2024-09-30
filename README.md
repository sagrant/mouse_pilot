# mouse_pilot
Bioinformatics workflow for gut microbiome analysis following a Proteomic Stable Isotope Probing experiment. This workflow is designed to process results from an experiment that investigated the effect of prebiotics on the mouse gut microbiome. 

Prebiotics are carbohydrates that are exclusively degraded by beneficial bacteria in the gut microbiome, and they improve host health by providing a nutrient source for these beneficial populations. We do not have a way to identify the specific microbial species that degrade prebiotics, which restricts our ability to understand how we could support those populations. To address this knowledge gap, we use a method called Proteomic Stable Isotope probing, where a prebiotic composed of stable isotope-labeled carbons is administered to live mice. That stable isotope will be assimilated by taxa that are actively degrading the prebiotic, which ultimately allows us to link the prebiotic substrate to the taxa that use it. This is done using a combination of metagenomics and proteomics. 

There are two sides of this analysis: 
1. Generation of metagenome annotations
3. Analysis of proteomics results

Overview of entire workflow:
![github_figure1](https://github.com/user-attachments/assets/fed510af-562b-40a2-8b10-82fdf71e015e)



A diagram of the metagenome annotations workflow is displayed below:
![github_figure2](https://github.com/user-attachments/assets/91e0f049-ef71-49d2-b973-faa36fe2c526)
The goal of this analysis is to assign taxonomic annotations to genes, contigs, and bins. Multiple softwares are used to increase our confidence that each annotation is correct. Our approach combines the Least Common Ancestor (LCA) and majority rule methods to assign annotations to genes, contigs, or bins.


Description of each step in the metagenome annotations workflow:
1. (a) Generate contig-level annotations with mmseqs software (b) Generate gene-level annotations with DIAMOND software.
2. Generate dummy .taxid file for gb_taxonomy_tools. Normal .taxid files have 3 columns: genbank ID, taxon ID, and count. A dummy .taxid simulates the first and third columns with irrelevant integers. The second column contains valid taxon IDs. This allows the researcher to retrieve all taxonomic ranks for each DIAMOND annotation without carrying out the entire gb_taxonomy_tools workflow.
3. Run gb_taxonomy_tools with dummy .taxid file.
4. Merge DIAMOND and gb_taxnonomy_tools output. The output file contains 4 columns: Gene ID, Protein ID, NCBI Taxon ID, Taxonomic Ranks.
5. Generate consensus annotations for each gene. 
6. Generate consensus annotations for each contig.
7. Combine all annotation data to facilitate comparison. Determine if annotations output by different softwares agree with one another.
8. Generate consensus annotations for each bin. Consensus annotations are provided at the gene level if there is no consensus within a bin.



