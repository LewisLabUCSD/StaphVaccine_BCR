# StaphVaccine_BCR
This study is for analyzing BCR data of staph vaccine (Cell Host and Microbe).

I. Script for generate Figure 4A
   compareClonotypes.R: this script is used to compare clone types and generates figure 4a and EXCEL file (Clonotypes & Sample).
   – Input file: filtered_contig_annotations_SampleName.csv (generated from Cell Ranger "vdj" output). (Note: SampleName = "NT_IsdB","LAC_alum","LAC_lsdB")
   – Output files: 
           (1) "compareClonotypes_aa.lsdB2.N=40.pdf" in the folder "results".
           (2) "compareClonotypes_aa.lsdB2.N=40.xlsx" in the folder "results".

II.Script for generate Figure 4B
   clononetwork.py: generates figure 4b. 
      – Input files: 
           (1) filtered_contig_annotations_SampleName.csv (generated from Cell Ranger "vdj" output. (Note: SampleName = "NT_IsdB","LAC_alum","LAC_lsdB")
           (2) the filtered_feature_bc_matrix.h5 (generated from Cell Ranger "count" output). 
      – Output file: "10_clonotype_network_5.png" in the folder "results".
