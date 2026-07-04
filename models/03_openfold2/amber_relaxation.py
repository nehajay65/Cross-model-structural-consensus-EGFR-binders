# From: https://github.com/deepmind/alphafold/issues/721

import json
import os
import pathlib
import pickle
import random
import sys
import time
import argparse
from typing import Dict

from absl import app
from absl import flags
from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model
from alphafold.relax import relax
import numpy as np

parser = argparse.ArgumentParser(description='Run Amber Relax (AlphaFold2 settings) on any structure')
parser.add_argument('-output_dir', type = str, help = 'Output path for the relaxed model')
parser.add_argument("input_dir", help="Directory containing unrelaxed .pdb files")

# An option to tell where the AF2 code is. 
# Also change the **default** to wherever your cloned AF2 is.
parser.add_argument("--af2_dir", default="home/postyr/alphafold", help="AlphaFold code directory")
 
args = parser.parse_args()
# This line below should put the AF code in your environment path
sys.path.append(args.af2_dir) 

def relax_with_amber(model_path,output_dir):
        RELAX_MAX_ITERATIONS = 0
        RELAX_ENERGY_TOLERANCE = 2.39
        RELAX_STIFFNESS = 10.0
        RELAX_EXCLUDE_RESIDUES = []
        RELAX_MAX_OUTER_ITERATIONS = 20
        amber_relax = relax.AmberRelaxation(use_gpu=True, max_iterations=RELAX_MAX_ITERATIONS, tolerance=RELAX_ENERGY_TOLERANCE, stiffness=RELAX_STIFFNESS, exclude_residues=RELAX_EXCLUDE_RESIDUES, max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS)
        model_path = pathlib.Path(model_path)
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Relaxing: {model_path}")

        with open(model_path, "r") as f:
            test_prot = protein.from_pdb_string(f.read())
                
        pdb_min, _, _ = amber_relax.process(prot=test_prot)
        
        if "unrelaxed" in model_path.name:
            output_name = model_path.name.replace("unrelaxed", "relaxed")
        else:
            output_name = model_path.stem + "_relaxed.pdb"
        
        output_path = output_dir / output_name

        with open(output_path, "w") as rel_f:
            rel_f.write(pdb_min)
        
        print(f"Saved: {output_path}")

def main():
    input_dir = pathlib.Path(args.input_dir)
    output_dir = pathlib.Path(args.output_dir)

    models = sorted(input_dir.glob("*.pdb"))

    if not models:
        print(f"No .pdb files found in {input_dir}")
    
    for model_path in models:
        # Skip if already relaxed
        if "unrelaxed" in model_path.name:
            output_name = model_path.name.replace("unrelaxed", "relaxed")
        else:
            output_name = model_path.stem + "_relaxed.pdb"
        
        if (output_dir / output_name).exists():
            print(f"Already done, skipping: {model_path.name}")
            continue

        try:
            relax_with_amber(model_path, output_dir)
        except Exception as e:
            print(f"FAILED: {model_path.name} — {e}")
            continue

if __name__ == "__main__":
    main()
