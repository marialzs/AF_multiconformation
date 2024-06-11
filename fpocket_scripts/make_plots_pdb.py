"""
Created on Mon Jun 12 17:07:33 2023

@author: masha.lzs
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Parse the input file names using argparse
parser = argparse.ArgumentParser(description='Plot distribution of pocket volumes.')
parser.add_argument('-f', help='Input file name for data', required=True)
parser.add_argument('-l', help='File name for list of bound structures', required=True)
args = parser.parse_args()

# Read the data file into a DataFrame
df = pd.read_csv(args.f, delim_whitespace=True, header=0)
print(df)
# Set the column names, adjust this based on the actual format of your file
df = df[['pdb', 'pock_vol']] 

# Read the list of bound structures
with open(args.l, 'r') as file:
    bound_structures = [line.strip() for line in file]
print(bound_structures)
bound_structures = bound_structures[1:]

# Create a new column in the DataFrame to indicate whether the structure is bound
df['is_bound'] = [1 if x.split('/')[1].split('.pdb')[0] in bound_structures else 0  for x in df['pdb']]
print(df)
# Filter the DataFrame for bound and unbound structures
df['is_bound'] = df['is_bound'].astype(bool)
bound_df = df[df['is_bound']]
unbound_df = df[~df['is_bound']]

# Plot the distribution of pock_vol values for bound and unbound structures
plt.figure(figsize=(8, 6))
plt.hist(unbound_df['pock_vol'].values, bins=20,  color='blue', edgecolor='black', alpha=0.5, label='Unbound')
plt.hist(bound_df['pock_vol'].values, bins=20, color='red', edgecolor='black', alpha=0.5, label='Bound')
plt.xlabel('Pocket Volume', fontsize=20)
plt.ylabel('Frequency', fontsize=20)
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 2200)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.legend(fontsize=14)
plt.savefig('pock_vol_distribution_overlay.png')
plt.show()
