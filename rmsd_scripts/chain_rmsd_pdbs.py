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


import Bio.PDB
from Bio.pairwise2 import align
import pandas as pd
from pymol import cmd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from matplotlib.colors import LinearSegmentedColormap
import math

#parse path
parser = argparse.ArgumentParser()
parser.add_argument('dir', help = 'directory with pdbs to look at')
parser.add_argument('schains')
parser.add_argument('close')
parser.add_argument('open')
parser.add_argument('path')
parser.add_argument('site', default = '0')
parser.add_argument('bound')
args = parser.parse_args()

# Function to save list to a text file
def save_list_to_file(data, filename):
    with open(filename, 'w') as file:
        for item in data:
            file.write(item + '\n')

#%%import csv with structure names

################################################################
###FETCH, ALIGN AND SAVE AS PDB ALL THE STRUCTURES IN THE LIST##
################################################################

#create a list of pdbs
pdb_list = [filename for filename in os.listdir(args.dir) if filename.endswith('.pdb')]


resis = args.schains.split(',')
atoms_to_be_aligned = [int(i) for i in resis]
#print(atoms_to_be_aligned)
exclude_selection_ref = f' and not resi '
for res in atoms_to_be_aligned:
    exclude_selection_ref = exclude_selection_ref + str(res) + '+'
exclude_selection_ref = exclude_selection_ref[:-1]


exclude_selection_set = f' and not resi '
for res in atoms_to_be_aligned:
    exclude_selection_set = exclude_selection_set + str(res+1) + '+'
exclude_selection_set = exclude_selection_set[:-1]


#fetch and align open one by one + save the structures
cmd.fetch(args.open, name = 'open')
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
open_chain = pdb_list[0].split('.')[-1]
open_model  = open_structure[0]

#read the PDB for open
close_structure = pdb_parser.get_structure("close", 'close/' + close_name)
close_chain = pdb_list[0].split('.')[-1]
close_model  = close_structure[0]



open_atoms = []
#iterate to append correct residues
for open_chain in open_model:   #put the chain you want to align
  # Iterate of all residues in each model in order to find proper atoms
  for open_res in open_chain:
    # Check if residue number ( .get_id() ) is in the list
    if open_res.get_id()[1] in atoms_to_be_aligned:
        sorted_open_atoms = sorted(open_res, key=lambda atom: atom.get_name())
        print('open', len(sorted_open_atoms))
        open_atoms.extend(sorted_open_atoms)


close_atoms = []
#iterate to append correct residues
for close_chain in close_model:   #put the chain you want to align
  # Iterate of all residues in each model in order to find proper atoms
  for close_res in close_chain:
    # Check if residue number ( .get_id() ) is in the list
    if close_res.get_id()[1] in (atoms_to_be_aligned):
        sorted_close_atoms = sorted(close_res, key=lambda atom: atom.get_name())
        print('closed', len(sorted_close_atoms))
        close_atoms.extend(sorted_close_atoms)

      
distance_to_closed = []
distance_to_open = []
ind = []

new_pdb_list = []


#run through all structures
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
                    print(len(sorted_sample_atoms))
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
                     new_pdb_list.append(structure)
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



#import bound txt
bound_info = pd.read_csv(args.bound)
bound_info = bound_info['Structure'].tolist()

#create logical list of shapes
shape = [1 if pdb[:-4] in bound_info else 0 for pdb in new_pdb_list] #for bound shape 1
#save scores in order
scores = pd.read_csv(args.path + '/structure_bs_strength.csv')
colormap = [scores[scores['Structure'] == pdb[19:-4]]['Binding Site 00' + args.site].iloc[0] for pdb in new_pdb_list]


#make figure
fig, ax = plt.subplots()
ax.set_box_aspect(1)
#Choose a colormap
cmap = plt.get_cmap('magma')

markers = ['o' if s == 1 else '^' for s in shape]

# Create a ScalarMappable for the colormap
cmap = LinearSegmentedColormap.from_list('blue_to_red', ['blue', 'red'])
norm = plt.Normalize(vmin=min(colormap), vmax=np.max(colormap))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Plot each point individually
for i in range(len(distance_to_closed)):
    plt.scatter(distance_to_closed[i], distance_to_open[i],edgecolors='black', s =50,linewidths=0.5, color=sm.to_rgba(colormap[i]), marker=markers[i])

# Add a colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=np.max(colormap)))
sm.set_array([])
cbar = plt.colorbar(sm,ticks = np.arange(0,np.max(colormap) + 1,5))
cbar.ax.tick_params(labelsize=14)
#Add labels and show the plot
plt.xlabel('Side Chain RMSD to Closed Structure (Å)', fontsize = 15)
plt.ylabel('Side Chain RMSD to Open Structure (Å)', fontsize = 15)

# Add a legend for shapes
bound_marker = plt.scatter([], [], marker='o', color='black', label='Bound')
unbound_marker = plt.scatter([], [], marker='^', color='black', label='Unbound')
plt.legend(handles=[bound_marker, unbound_marker], loc='upper left', fontsize = 12)
plt.tight_layout()

#ax.set_xlim(-0.5,4)
#ax.set_ylim(-0.5,4)
#ax.set_xticks(np.arange(0,5,1))
#ax.set_yticks(np.arange(0,5,1))
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)

plt.show()
r_c, pi_c = stats.pearsonr(distance_to_closed, colormap)
r_o, pi_o = stats.pearsonr(distance_to_open, colormap)

print(f'The correlation score of RMSD to closed structure and FTMove scores is: {r_c:.3f}')
print(f'The correlation score of RMSD to open structure and FTMove scores is: {r_o:.3f}')



# Initialize variables to track the indices and values
closed_str = []
open_str = []
other = []
# Iterate over the index and value pairs from both x and y
for index, (xi, yi) in enumerate(zip(distance_to_closed, distance_to_open)):
    distance_x = math.sqrt((1.95 - xi) ** 2 + (0 - yi) ** 2)
    distance_y = math.sqrt((0 - xi) ** 2 + (1.95 - yi) ** 2)
    if distance_x < 2.76/2:
        open_str.append(new_pdb_list[index])

    elif distance_y < 2.76/2:
        closed_str.append(new_pdb_list[index])
    else:
        other.append(new_pdb_list[index])

print('Closed Structures', len(closed_str))
print('Open Structures', len(open_str))
print('Other', len(other))

# Save each list to a corresponding text file
save_list_to_file(closed_str, 'closed_structures.txt')
save_list_to_file(open_str, 'open_structures.txt')
save_list_to_file(other, 'other_structures.txt')

