
import os
import os.path
import re
import hashlib
import random
import sys
import string
import urllib.request
import subprocess
import fileinput
import argparse
import warnings
import logging
from Bio import BiopythonDeprecationWarning
from pathlib import Path
import matplotlib.pyplot as plt
from colabfold.download import download_alphafold_params, default_data_dir
from colabfold.utils import setup_logging
from colabfold.batch import get_queries, run, set_model_type
from colabfold.plot import plot_msa_v2
from colabfold.colabfold import plot_protein
import numpy as np
import py3Dmol
import glob
from colabfold.colabfold import plot_plddt_legend
from colabfold.colabfold import pymol_color_list, alphabet_list
from IPython.display import display, HTML
import base64
from html import escape
import torch
import shutil
from PIL import Image
import logging
from Bio.PDB import PDBParser, PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from daqrefine import Daqrefine



def search_files(directory, extension):
    return [os.path.join(root, file)
            for root, dirs, files in os.walk(directory)
            for file in files if file.endswith(extension)]

def get_arguments():
    parser = argparse.ArgumentParser(description='STEP-1: Input Protein Sequence and DAQ result file')

    # Add arguments
    parser.add_argument('--str_mode', type=str, default='strategy 2',
                        help='Select the DAQ-refine strategy. Choices are Vanilla AF2, strategy 1, and strategy 2.')
    
    # MENSMMFISRSLRRPVTALNCNLQSVRTVIYLHKGPRINGLRRDPESYLRNPSGVLFTEVNAKECQDKVRSILQLPKYGINLSNELILQCLTHKSFAHGSKPYNEKLNLLGAQFLKLQTCIHSLKNGSPAESCENGQLSLQFSNLGTKFAKELTSKNTACTFVKLHNLGPFIFWKMRDPIKDGHINGETTIFASVLNAFIGAILSTNGSEKAAKFIQGSLLDKEDLHSLVNIANENVASAKAKISDKENKAFL
    parser.add_argument('--query_sequence', type=str, default='',
                        help='Input target protein sequence.',required=True)

    parser.add_argument('--jobname', type=str, default='',
                        help='Job name for the task.',required=True)
    
    parser.add_argument('--num_relax', type=int, default=0,
                        help='Number of models to use.')
    
    parser.add_argument('--template_mode', type=str, default='none',
                        help='Template mode: none, pdb70, or custom.')
    
    parser.add_argument('--pdb_input_path', type=str, default='',
                        help='Path to the DAQ-score output file (if applicable).')
    
    parser.add_argument('--cust_msa_path', type=str, default='none',
                        help='Path to the custom MSA file (if applicable).')

    parser.add_argument('--input_path', type=str, default='',
                        help='Path to the directory of the input files.',required=True)

    parser.add_argument('--output_path', type=str, default='',
                        help='Path to the directory of the output files.',required=True)

    parser.add_argument('--resolution', type=str, default='3.43',
                        help='Specify the resolution in the relaxation',required=True)

    args = parser.parse_args()

    return args

def main():
    # Get arguments (this function needs to be implemented based on the original code)
    args = get_arguments()

    # run s1
    print("INFO: STEP-1: Running strategy 1")
    args.str_mode = 'strategy 1'

    s1 = Daqrefine(
        args=args
    )
    # clean up the directory
    s1.clean_up(args)
    # run the s1 modeling
    s1.run_modeling()

    # run s2
    # run vanilla alphafold to get msa file
    print("INFO: STEP-2: Running Vanilla AF2")
    vanilla_af2_result = None
    args.str_mode = 'Vanilla AF2'
    vanilla_af2_result = Daqrefine(
        args=args
    )
    # run vanilla alphafold modeling
    vanilla_af2_result.run_modeling()

    print("STEP-3: Running strategy 2")
    args.str_mode = 'strategy 2'
    a3m_files = search_files(vanilla_af2_result.result_dir, '.a3m')
    args.cust_msa_path = a3m_files[0]

    s2 = Daqrefine(
        args=args
    )
    
    # run the s2 modeling process
    s2.run_modeling()


if __name__ == '__main__':
    main()

