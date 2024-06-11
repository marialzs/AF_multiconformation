#!/bin/bash

"""
Created on Mon Jun 12 17:07:33 2023

@author: masha.lzs
"""

ref_pdb="$1"
ref_lig="$2"
struct_dir="$3"

#activate environment w pymol api
conda activate pymol_env

echo "Creating directory: ${struct_dir%/}_w_lig"

mkdir "${struct_dir%/}_w_lig"
#align all structs w ref and add ligand
python ../align_and_addlig.py "$ref_pdb" "$ref_lig" "$struct_dir" "${struct_dir%/}_w_lig"


#make list of structure vs ligand (same ligand in this case)
python ../make_ligand_list.py "${struct_dir%/}_w_lig" "${ref_lig%.pdb}"

conda deactivate
conda activate fpocket_env

#run dpocket
dpocket -f ligand_pairs.txt

#plot and save the results
python ../plot_scores.py

#move results to a folder
mkdir "${struct_dir%/}_dpocket"
mv dpout* "${struct_dir%/}_dpocket/"
mv ligand_pairs.txt "${struct_dir%/}_dpocket/"

echo "Analysis done"
