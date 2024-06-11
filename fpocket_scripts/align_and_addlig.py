"""
Created on Mon Jun 12 17:07:33 2023

@author: masha.lzs
"""

import os
import pymol
from pymol import cmd
import argparse

def align_and_add_ligand(ref_pdb, ligand_pdb, struct_dir, output_dir):
    """
    Aligns PDB structures to a reference structure and adds a ligand to each aligned structure.

    Args:
    ref_pdb (str): Path to the reference PDB file.
    ligand_pdb (str): Path to the ligand PDB file.
    struct_dir (str): Directory containing the PDB structures to be aligned.
    output_dir (str): Directory where the modified PDB files will be saved.
    """

    #Load the reference structure and ligand
    cmd.load(ref_pdb, 'ref_structure')
    cmd.load(ligand_pdb, 'ligand')

    #Iterate over the PDB files in the structure directory
    for pdb_file in os.listdir(struct_dir):
        if pdb_file.endswith('.pdb'):
            struct_path = os.path.join(struct_dir, pdb_file)
            cmd.load(struct_path, 'target_structure')

            #Align the target structure to the reference structure
            cmd.align('target_structure', 'ref_structure')

            #Remove all heteroatoms from the target structure
            cmd.remove('target_structure and het')

            #Add the ligand to the target structure
            cmd.create('modified_structure', 'target_structure + ligand')

            #Save the modified structure to the output directory
            output_path = os.path.join(output_dir, pdb_file)
            cmd.save(output_path, 'modified_structure')

            #Delete the structures from the PyMOL session
            cmd.delete('target_structure')
            cmd.delete('modified_structure')


#parse path
parser = argparse.ArgumentParser()
parser.add_argument('ref_pdb')
parser.add_argument('ligand_pdb')
parser.add_argument('struct_dir')
parser.add_argument('output_dir')
args = parser.parse_args()


align_and_add_ligand(args.ref_pdb, args.ligand_pdb, args.struct_dir, args.output_dir)

