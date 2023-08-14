# SpecLoop
SpecLoop is a method for predicting cell type-specific chromatin loops via transcription factor cooperation.
# Motivation
Cell-type-Specific Chromatin Loops (CSCLs) coincide with different chromatin states, downstream gene regulation, and cell type-specific functions. However, the establishment mechanisms of CSCL remains unclear and experimentally measuring genome-wide CSCL is technically challenging, costly, and time-consuming. We notice that genome-wide transcriptomic and epigenomic data for various tissues and cell types are rapidly accumulated and publicly accessible in many databases, such as the Encyclopedia of DNA Elements (ENCODE) and Roadmap Epigenomic database. Therefore, developing computational approach to predict CSCL for the cell types of interest from the widely available multi-omics data is desirable.
# Method
We develop SpecLoop, a network regularization-based machine learning framework, to explore the role of transcription factors (TFs) cooperation in the formation of CSCL. SpecLoop predicts CSCLs and identifies TF cooperations by integrating widely available gene expression, chromatin accessibility, sequence, protein-protein interaction, and TF binding motif data. We take advantage of network regularization to choose TF and TF complex related features weight efficiently by reconstructing TF protein-protein interaction and co-expression network in an unsupervised manner. SpecLoop accurately predicted CSCL in a range of cell types, and outperformed other sparsity regularized methods in both prediction accuracy and feature selection. Moreover, SpecLoop identified TFs and TF combinations as potential participants in the formation of CSCL chromatin loops. This provides valuable biological insight into understanding the molecular mechanisms underlying the establishment of CSCLs.
# SpecLoop Source code
SpecLoop
Version 1.0 Last updated: August 8, 2023
# Requirements
R 4.0.4
R package: glmnet, Matrix, foreach, pROC, and PRROC.
# Reference
Lixin Ren, Wanbiao Ma, and Yong Wang, SpecLoop predicts cell type-specific chromatin loop via transcription factor cooperation. In submission.
