# mouse_pilot
Bioinformatics workflow for mouse gut microbiome analysis

Premise:
This workflow is designed to process results from a gut microbiome experiment that looked at the effect of prebiotics on the mouse gut microbiome. 

Prebiotics are carbohydrates that are exclusively degraded by beneficial bacteria in the gut microbiome, and they improve host health by providing a nutrient source for these beneficial populations. We do not have a way to detect which microbial species actively degrade prebiotics, which restricts our ability to understand how we could better support those populations with these compounds. For this, we use a method called Proteomic Stable Isotope probing, where a prebiotic composed of stable isotope-labeled carbons is administered to live mice. That stable isotope will be taken up by taxa that are actively degrading the prebiotic, which ultimately allows us to link the prebiotic substrate to the taxa that assimilate it. This is done using a combination of metagenomics and proteomics. 

There are two sides of this analysis: 
1. Generation of metagenome annotations
2. Analysis of proteomics results

Overview of entire workflow:
![github_figure1](https://github.com/user-attachments/assets/fed510af-562b-40a2-8b10-82fdf71e015e)



A diagram of the metagenome annotations workflow is displayed below:
![github_figure2](https://github.com/user-attachments/assets/91e0f049-ef71-49d2-b973-faa36fe2c526)
Description of each step in the metagenome annotations workflow:
1a. Generate contig-level annotations with mmseqs software
1b. Generate gene-level annotations with DIAMOND software
2. Generate dummy .taxid file for gb_taxonomy_tools. Normal .taxid files have 3 columns: genbank ID, taxon ID, and count. A dummy .taxid simulates the first and third columns with superfluous integers. The second column contains valid taxon IDs. This allows the researcher to get all taxonomic ranks for each DIAMOND annotation without the 



