#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:17:56 2023

@author: masha.lzs

This code takes in an rmsd matrix and clusters indexes to give out 
representatives for each
"""

import argparse
import numpy as  np
import math


parser = argparse.ArgumentParser(description='Process a file with a given radius.')

# Add arguments
parser.add_argument('-f', '--file', type=str, help='File path')
parser.add_argument('-r', '--radius', type=float, help='Radius')
parser.add_argument('-t', '--type', type=str, help='reps or all')

# Parse the command-line arguments
args = parser.parse_args()

mat = np.loadtxt(args.file)

try:
    mat.shape[1]    
except IndexError:
    newmat = mat.reshape([int(math.sqrt(len(mat))), int(math.sqrt(len(mat)))], order = 'c')
rad = args.radius
#%%
clusters = {}
visited = []
#%%
dicind = {}  
count = 0

while len(dicind) > 0 or len(visited) == 0:
    dicind = {}  
    diclen = {} 
    for i,col in enumerate(np.arange(0, len(newmat), 1)):
        indexes = []
        if col in visited:
            pass
        else:
            for j,row in enumerate(np.arange(0, len(newmat), 1)):
                if newmat[row,col] <= rad and row not in visited:
                    indexes.extend([j])
            dicind.update({col: indexes})
            diclen.update({col: len(indexes)})

    #find index for highest cluster
    if len(diclen) > 0:
        biggest_rep = max(diclen, key = lambda k: diclen[k]) 
        clusters.update({biggest_rep: dicind.get(biggest_rep)})
        visited.extend(dicind.get(biggest_rep))
        if biggest_rep not in visited:
            visited.extend([biggest_rep]) 

if args.type == 'reps': 
    for i,rep in enumerate(clusters.keys()):
        print(rep)   
else:
    print(clusters.values())   
    


 
    

            
    
    
            






    
    
