# StaphVaccine_BCR
This study is for analyzing BCR data of staph vaccine (Cell Host and Microbe).

## I. Script for generate Figure 4A
1. compareClonotypes.R: this script is used to compare clone types and generates figure 4a and EXCEL file (Clonotypes & Sample).
   - Input file: filtered_contig_annotations_SampleName.csv (generated from Cell Ranger "vdj" output). 
     - Note: SampleName = "NT_IsdB","LAC_alum","LAC_lsdB"
   - Output files: 
     - "compareClonotypes_aa.lsdB2.N=40.pdf" in the folder "results".
     - "compareClonotypes_aa.lsdB2.N=40.xlsx" in the folder "results".

## II. Script for generate Figure 4B
2. clononetwork.py: generates figure 4b. 
   - Input files: 
     - filtered_contig_annotations_SampleName.csv (generated from Cell Ranger "vdj" output). 
     - the filtered_feature_bc_matrix.h5 (generated from Cell Ranger "count" output). 
   - Output file: "10_clonotype_network_5.png" in the folder "results".
