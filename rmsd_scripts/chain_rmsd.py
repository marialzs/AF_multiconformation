#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 17:40:21 2023

@author: masha.lzs
"""

# -*- coding: utf-8 -*-
# The MIT License
# 
# Copyright (c) 2010-2016 Anders S. Christensen
# 
#code modified from the owner above
#some parts of the code taken from license above

import matplotlib.lines as mlines
import Bio.PDB
from Bio.pairwise2 import align
import pandas as pd
from pymol import cmd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import math

#parse path
parser = argparse.ArgumentParser()
parser.add_argument('dir', help = 'directory with pdbs to look at')
parser.add_argument('schains')
parser.add_argument('close')
parser.add_argument('open')
parser.add_argument('analysis', default = 'single')
parser.add_argument('site', default = 'none')
args = parser.parse_args()


#%%import csv with structure names

################################################################
###FETCH, ALIGN AND SAVE AS PDB ALL THE STRUCTURES IN THE LIST##
################################################################
#find rasidue range
resis = args.schains.split(',')
atoms_to_be_aligned = [int(i) for i in resis]
exclude_selection_ref = f' and not resi '
for res in atoms_to_be_aligned:
    exclude_selection_ref = exclude_selection_ref + str(res) + '+'
exclude_selection_ref = exclude_selection_ref[:-1]


exclude_selection_set = f' and not resi '
for res in atoms_to_be_aligned:
    exclude_selection_set = exclude_selection_set + str(res) + '+'
exclude_selection_set = exclude_selection_set[:-1]

#create a list of pdbs
pdb_list = [filename for filename in os.listdir(args.dir) if filename.endswith('.pdb')]
if args.analysis != 'single':
    analysis_list = [filename for filename in os.listdir(args.analysis) if filename.endswith('.pdb')]

#fetch and align open one by one + save the structures
cmd.fetch(args.open,name = 'open')
open_name = args.open + '_al.pdb'

#loop through structures for open
for structure in pdb_list:
    name = structure + '_al.pdb'
    cmd.load(args.dir + structure, structure)
    cmd.align(structure+exclude_selection_set, 'open'+exclude_selection_ref)
    
    cmd.save('open/' + name, selection = structure)

cmd.save('open/' + open_name, selection = 'open' )

#delete session
cmd.delete('*')
    
#fetch and align closed one by one + save the structures
cmd.fetch(args.close, name = 'close')
close_name = args.close + '_al.pdb'

#loop through structures for  closed
for structure in pdb_list:
    
    name = structure + '_al.pdb'
    cmd.load(args.dir + structure, structure)
    cmd.align(structure+exclude_selection_set, 'close'+exclude_selection_ref)
    
    cmd.save('close/' + name, selection = structure)

cmd.save('close/' + close_name, selection = 'close' )


#read the PDB for open
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
open_structure = pdb_parser.get_structure("open", 'open/' + open_name)
open_chain = args.open.split('.')[-1]
open_model  = open_structure[0]

#read the PDB for open
close_structure = pdb_parser.get_structure("close", 'close/' + close_name)
close_chain = args.close.split('.')[-1]
close_model  = close_structure[0]



open_atoms = []
#iterate to append correct residues
for open_chain in open_model:   #put the chain you want to align
  # Iterate of all residues in each model in order to find proper atoms
  for open_res in open_chain:
    print(open_res.get_id())
    # Check if residue number ( .get_id() ) is in the list
    if open_res.get_id()[1] in atoms_to_be_aligned:
        sorted_open_atoms = sorted(open_res, key=lambda atom: atom.get_name())
        open_atoms.extend(sorted_open_atoms)

close_atoms = []
#iterate to append correct residues
for close_chain in close_model:   #put the chain you want to align
  # Iterate of all residues in each model in order to find proper atoms
  for close_res in close_chain:
    # Check if residue number ( .get_id() ) is in the list
    if close_res.get_id()[1] in (atoms_to_be_aligned):
        sorted_close_atoms = sorted(close_res, key=lambda atom: atom.get_name())
        close_atoms.extend(sorted_close_atoms)
      
#find residue range
atoms_to_be_aligned = [int(i) for i in resis]


#find residue range
distance_to_closed = []
distance_to_open = []
centr_distance_open = []
centr_distance_closed = []

for path in ['open', 'close']:      
    for structure in pdb_list:
        
       	try:
            sample_atoms = []
            sample_structure = pdb_parser.get_structure("sample", path + '/' + structure + '_al.pdb')
            sample_chain = structure.split('.')[-1]
            sample_model = sample_structure[0]
            
            # Iterate of all chains in the model in order to find all residues
            for sample_chain in sample_model:
              for sample_res in sample_chain:
                if sample_res.get_id()[1] in atoms_to_be_aligned:
                    sorted_sample_atoms = sorted(sample_res, key=lambda atom: atom.get_name())
                    sample_atoms.extend(sorted_sample_atoms)

            if len(sample_atoms) == len(open_atoms):
                
                 errorsum = 0
                 
                 if path == 'open':
                     for i in range(len(sample_atoms)):
                         for j in range(3):
                             try:
                                 
                                 errorsum += (open_atoms[i].get_vector()[j] - sample_atoms[i].get_vector()[j])**2
                             except IndexError:
                                 pass
                     RMSD = (errorsum/len(sample_atoms))**(1/2)
                     distance_to_open.append(RMSD)
                 else:
                      for i in range(len(sample_atoms)):
                         for j in range(3):
                             try:
                                 errorsum += (close_atoms[i].get_vector()[j] - sample_atoms[i].get_vector()[j])**2
                             except IndexError:
                                 pass
                      RMSD = (errorsum/len(sample_atoms))**(1/2)
                      distance_to_closed.append(RMSD)
        except KeyError:
            pass

    if args.analysis != 'single':
        for structure in analysis_list:
        
            try:
                print(structure)
                sample_atoms = []
                sample_structure = pdb_parser.get_structure("sample", path + '/' + structure + '_al.pdb')
                sample_chain = structure.split('.')[-1]
                sample_model = sample_structure[0]

            # Iterate of all chains in the model in order to find all residues
                for sample_chain in sample_model:
                  for sample_res in sample_chain:
                    if sample_res.get_id()[1] in atoms_to_be_aligned:
                        sorted_sample_atoms = sorted(sample_res, key=lambda atom: atom.get_name())
                        sample_atoms.extend(sorted_sample_atoms)

                if len(sample_atoms) > 0:

                     errorsum = 0

                     if path == 'open':
                         for i in range(len(sample_atoms)):
                             for j in range(3):
                                 errorsum += (open_atoms[i].get_vector()[j] - sample_atoms[i].get_vector()[j])**2
                         RMSD = (errorsum/len(sample_atoms))**(1/2)
                         centr_distance_open.append(RMSD)
                     else:
                          for i in range(len(sample_atoms)):
                             for j in range(3):
                                 errorsum += (close_atoms[i].get_vector()[j] - sample_atoms[i].get_vector()[j])**2
                          RMSD = (errorsum/len(sample_atoms))**(1/2)
                          centr_distance_closed.append(RMSD)

            except KeyError:
                pass


#check outliers
for i,item in enumerate(distance_to_closed):
    for j,item2 in enumerate(distance_to_open):
        if i == j and item > 1 and item < 2 and item2 > 2 and item2 < 2.5:
            print(pdb_list[i], item, item2)
#plot a scatterplot
fig, ax = plt.subplots()

ax.set_box_aspect(1)
scatter = ax.scatter(distance_to_closed,distance_to_open,s = 50, color = 'grey',
               alpha=0.5)

if args.analysis != 'single':
    #import csv
    scores_in_order = []
    if args.site != 'none':
        csv_files = [file for file in os.listdir(args.analysis) if file.endswith('.csv')]
        binding_sites = pd.read_csv(args.analysis + '/' + csv_files[0])
        bs_scores = binding_sites['Binding Site 00' + args.site]
        for i,structure in enumerate(analysis_list):
         score = int(binding_sites['Binding Site 00' + args.site][binding_sites['Structure'] == structure.split('.')[0]])
         scores_in_order.append(score)
    else:
        score = 30
    # Define a custom colormap from blue to red
    cmap = LinearSegmentedColormap.from_list('blue_to_red', ['blue', 'red'])

# Normalize the colormap based on the range of scores
    norm = plt.Normalize(vmin=0, vmax=20)

# Create a ScalarMappable object with the normalized colormap
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

# Create the scatter plot with the custom colormap
    sc = plt.scatter(centr_distance_closed, centr_distance_open, s=100, cmap=cmap, c=scores_in_order, norm=norm)

# Add a colorbar using the ScalarMappable object
    cbar = plt.colorbar(sm, ticks=np.arange(0, 20 + 1, 5))
    cbar.ax.tick_params(labelsize=14)
handle1 = mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=5, alpha=0.5, label='AF Ensembles')
#handle2 = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, markeredgecolor='black', markeredgewidth=1, label='Cluster Centroids')

# Adding the legend
fig.legend(handles=[handle1],fontsize=12, 
           loc = 'upper right', bbox_to_anchor=((0.8, 0.95)))


#plt.title(args.schains+ " Residues RMSD")
plt.xlabel('Side Chain RMSD to Closed Structure (Å)', fontsize = 15)
plt.ylabel('Side Chain RMSD to Open Structure (Å)', fontsize = 15)
plt.xlim(-0.5,11)
plt.ylim(-0.5,9)
plt.xticks(np.arange(0,12,1))
plt.yticks(np.arange(0,10,1))
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.tight_layout()

plt.show()            

# Initialize variables to track the indices and values
best_index_x = None
best_index_y = None
shortest_dist_x = 40
shortest_dist_y = 40

# Iterate over the index and value pairs from both x and y
for index, (xi, yi) in enumerate(zip(distance_to_closed, distance_to_open)):
    distance_x = math.sqrt((6.5 - xi) ** 2 + (0 - yi) ** 2)
    distance_y = math.sqrt((0 - xi) ** 2 + (6.5 - yi) ** 2)
    if distance_x < shortest_dist_x:
        shortest_dist_x = distance_x
        best_index_x = index

    if distance_y < shortest_dist_y:
        shortest_dist_y = distance_y
        best_index_y = index
    
# Output the results
if best_index_x is not None and best_index_y is not None:
    # Calculate the Euclidean distance between the two points
    distance = math.sqrt((distance_to_closed[best_index_x] - distance_to_closed[best_index_y]) ** 2 + (distance_to_open[best_index_x] - distance_to_open[best_index_y]) ** 2)
    print(f"The distance between these two points is: {distance}")
else:
    print("No single data points match the required criteria.")

