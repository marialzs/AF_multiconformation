# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 17:07:33 2023

@author: masha.lzs
"""

from pymol import cmd
import numpy as np

filenames = []

#create a list of filenames
for i in range(1,6):
    for j in range(0, 100):
        if j < 10:
            filenames.append("unrelaxed_model_" + str(i) + "_ptm_pred_" + str(j) + ".pdb")
        else:
            filenames.append("unrelaxed_model_" + str(i) + "_ptm_pred_" + str(j) + ".pdb")

#load structs
for i,name in enumerate(filenames):
    cmd.load("../ribo_a/AF_ensembles/"+ name,name[:-4]) #modify this line


#calculate rmsd and create matrix
rmsd_mat = np.zeros([len(filenames), len(filenames)])

for i,name in enumerate(filenames):
    row = np.zeros([1, len(filenames)])    
    for j, other in enumerate(filenames):
        if i == j:
            pass
        else:
            rmsd, *stuff = cmd.align(name[:-4],other[:-4])
            row[0,j] = rmsd 
    rmsd_mat[i,:] = row

rmsd_mat = rmsd_mat.flatten('F')

#write in doc
with open("rmsd_riboa.txt", "w+") as file: #modify this line
    for item in rmsd_mat:
        file.write(f"{item}\n")
        

        


