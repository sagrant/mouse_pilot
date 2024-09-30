# mouse_pilot
Bioinformatics workflow for mouse gut microbiome analysis

Premise:
This workflow is designed to process results from a gut microbiome experiment that looked at the effect of prebiotics on the mouse gut microbiome. 

Prebiotics are carbohydrates that are exclusively degraded by beneficial bacteria in the gut microbiome, and they improve host health by providing a nutrient source for these beneficial populations. We do not have a way to detect which microbial species actively degrade prebiotics, which restricts our ability to understand how we could better support those populations with these compounds. For this, we use a method called Proteomic Stable Isotope probing, where a prebiotic composed of stable isotope-labeled carbons is administered to live mice. That stable isotope will be taken up by taxa that are actively degrading the prebiotic, which ultimately allows us to link the prebiotic substrate to the taxa that assimilate it. This is done using a combination of metagenomics and proteomics. 

There are two sides of this analysis: 
1. Generation of metagenome annotations
2. Analysis of proteomics results

