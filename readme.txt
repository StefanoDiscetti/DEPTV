Sample code supporting the scientific paper titled: 
 - "From sparse data to high-resolution fields: ensemble particle modes as a basis for high-resolution flow characterization"
 - Authors: J. Cortina Fernández, C. Sanmiguel Vila, A. Ianiro, S. Discetti.
 - Corresponding author: sdiscett@ing.uc3m.es

This folder contains two different MATLAB codes:
 1: main_debug.m - Sample code comparing PIV and DEPTV algorithms with the true DNS flow field
 2: main.m       - Sample code comparing PIV and DEPTV algorithms, intended to be modified by the user

The codes use a k-d tree function from Andrea Tagliasacchi: https://github.com/ataiya/kdtree
The DNS dataset has been modified from the Dynamic Mode Decomposition book: http://dmdbook.com/
Some plots are generated using a subtightplot function from Felipe G. Nievinski