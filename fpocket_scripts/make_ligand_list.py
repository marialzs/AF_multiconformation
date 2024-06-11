"""
Created on Mon Jun 12 17:07:33 2023

@author: masha.lzs
"""


import os
import argparse

def create_pdb_ligand_list(dir_path, ligand_id):
    """
    Creates a text file listing the path to each PDB file in the directory and the ligand ID.

    Args:
    dir_path (str): Path to the directory containing PDB files.
    ligand_id (str): Ligand ID to be associated with each PDB file.
    output_file (str): Path to the output text file.
    """
    with open('ligand_pairs.txt', 'w') as f:
        for file_name in os.listdir(dir_path):
            if file_name.endswith('.pdb'):
                pdb_path = os.path.join(dir_path, file_name)
                f.write(f'{pdb_path} {ligand_id}\n')

parser = argparse.ArgumentParser()
parser.add_argument('dir_path')
parser.add_argument('ligand_id')
args = parser.parse_args()


create_pdb_ligand_list(args.dir_path, args.ligand_id)

