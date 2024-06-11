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
parser.add_argument('start')
parser.add_argument('stop')
parser.add_argument('close')
parser.add_argument('open')
parser.add_argument('analysis', default = 'single')
parser.add_argument('site', default = 'none')
args = parser.parse_args()


#%%import csv with structure names

################################################################
###FETCH, ALIGN AND SAVE AS PDB ALL THE STRUCTURES IN THE LIST##
################################################################

exclude_range = (args.start, args.stop)

#create a list of pdbs
pdb_list = [filename for filename in os.listdir(args.dir) if filename.endswith('.pdb')]
if args.analysis != 'single':
    analysis_list = [filename for filename in os.listdir(args.analysis) if filename.endswith('.pdb')]

#fetch and align open one by one + save the structures
cmd.fetch(args.open, name = 'open')
open_name = args.open + '_al.pdb'

#loop through structures for open
for structure in pdb_list:
     
    name = structure + '_al.pdb'
    cmd.load(args.dir + structure, structure)
    exclude_selection = f' and not resi {exclude_range[0]}-{exclude_range[1]}'
    cmd.align(structure + exclude_selection, 'open'+ exclude_selection)
    
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
    exclude_selection = f' and not resi {exclude_range[0]}-{exclude_range[1]}'

    cmd.align(structure + exclude_selection, 'close'+ exclude_selection)
    
    cmd.save('close/' + name, selection = structure)

cmd.save('close/' + close_name, selection = 'close' )

#find rasidue range
start_id = int(args.start)
end_id = int(args.stop)
atoms_to_be_aligned = range(start_id, end_id + 1)

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
    # Check if residue number ( .get_id() ) is in the list
    if open_res.get_id()[1] in atoms_to_be_aligned:
      # Append CA atom to list
      open_atoms.append(open_res['CA'])
      
      
close_atoms = []
#iterate to append correct residues
for close_chain in close_model:   #put the chain you want to align
  # Iterate of all residues in each model in order to find proper atoms
  for close_res in close_chain:
    # Check if residue number ( .get_id() ) is in the list
    if close_res.get_id()[1] in atoms_to_be_aligned:
      # Append CA atom to list
      close_atoms.append(close_res['CA'])
      
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
                  sample_atoms.append(sample_res['CA'])
                  
            if len(sample_atoms) > 0:
                
                
                 errorsum = 0
                 if path == 'open':
                     for i in range(len(sample_atoms)):
                         for j in range(3):
                             errorsum += (open_atoms[i].get_vector()[j] - sample_atoms[i].get_vector()[j])**2
                     RMSD = (errorsum/len(sample_atoms))**(1/2)
                     distance_to_open.append(RMSD)
                 else:
                      for i in range(len(sample_atoms)):
                         for j in range(3):
                             errorsum += (close_atoms[i].get_vector()[j] - sample_atoms[i].get_vector()[j])**2
                      RMSD = (errorsum/len(sample_atoms))**(1/2)
                      distance_to_closed.append(RMSD)
            else:
                if path == 'open':
                    ind = pdb_list.index(structure)
                    color_map.pop(ind) 
        except KeyError:
            pass

    if args.analysis != 'single':
        for structure in analysis_list:
        
            try:
                sample_atoms = []
                sample_structure = pdb_parser.get_structure("sample", path + '/' + structure + '_al.pdb')
                sample_chain = structure.split('.')[-1]
                sample_model = sample_structure[0]

            # Iterate of all chains in the model in order to find all residues
                for sample_chain in sample_model:
                  for sample_res in sample_chain:
                    if sample_res.get_id()[1] in atoms_to_be_aligned:
                      sample_atoms.append(sample_res['CA'])

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
print(centr_distance_closed)
#check outliers
for i,item in enumerate(distance_to_open):
    if item > 5  and distance_to_closed[i] > 5:
        print(pdb_list[i],distance_to_closed[i], item)
                    
#plot a scatterplot
sp_names = ['Unbound', 'Bound']
#sp_names = ['AF Ensembles', 'Cluster Centroids']
fig, ax = plt.subplots()

ax.set_box_aspect(1)
scatter = ax.scatter(distance_to_closed,distance_to_open,s = 50,color = 'grey',
               alpha=0.5)
#c = color_map,cmap = 'coolwarm')
#fig.legend(handles=scatter.legend_elements()[0], 
#           labels=sp_names, loc = 'upper right', bbox_to_anchor=((0.8, 0.95)))
#plt.colorbar()
if args.analysis != 'single':
    #import csv
    scores_in_order = []
    if args.site != 'none':
        csv_files = [file for file in os.listdir(args.analysis) if file.endswith('.csv')]
        binding_sites = pd.read_csv(args.analysis + '/' + csv_files[0])
        bs_scores = binding_sites['Binding Site 00' + args.site]
        for i,structure in enumerate(analysis_list):
         score = int(binding_sites['Binding Site 00' + args.site][binding_sites['Structure'] == structure[19:-4]])
         scores_in_order.append(score)
    else:
        score = 30
    # Define a custom colormap from blue to red
    cmap = LinearSegmentedColormap.from_list('blue_to_red', ['blue', 'red'])

# Normalize the colormap based on the range of scores
    norm = plt.Normalize(vmin=0, vmax=21)

# Create a ScalarMappable object with the normalized colormap
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

# Create the scatter plot with the custom colormap
    sc = plt.scatter(centr_distance_closed, centr_distance_open, s=100, cmap=cmap, c=scores_in_order, norm=norm)

# Add a colorbar using the ScalarMappable object
    cbar = plt.colorbar(sm, ticks=np.arange(0, 21 + 1, 5))
    cbar.ax.tick_params(labelsize=14)
handle1 = mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=5, alpha=0.5, label='AF Ensembles')  # Set alpha for transparency
#handle2 = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, markeredgecolor='black', markeredgewidth=1, label='Cluster Centroids')  # Set marker edge color and width

# Adding the legend
fig.legend(handles=[handle1],fontsize=12, 
           loc = 'upper right',bbox_to_anchor=((0.8, 0.95)))

#plt.title(args.start + ' - ' + args.stop + ' alpha-carbon segment  RMSD')
ax.set_xlabel('Loop RMSD to Closed Structure (Å)', fontsize = 16)
ax.set_ylabel('Loop RMSD to Open Structure (Å)', fontsize = 16)
ax.set_xlim(-0.5,8)
ax.set_ylim(-0.5,7)
ax.set_xticks(np.arange(0,9,1))
ax.set_yticks(np.arange(0,8,1))
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
plt.tight_layout()

plt.show()            


# Initialize variables to track the indices and values
best_index_x = None
best_index_y = None
shortest_dist_x = 40
shortest_dist_y = 40

# Iterate over the index and value pairs from both x and y
for index, (xi, yi) in enumerate(zip(distance_to_closed, distance_to_open)):
    distance_x = math.sqrt((5.6 - xi) ** 2 + (0 - yi) ** 2)
    distance_y = math.sqrt((0 - xi) ** 2 + (5.6 - yi) ** 2)
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

