import argparse
import torch
import os

def argparser():
    parser = argparse.ArgumentParser()
    #for input file
    parser.add_argument("--log_folder_path", default="log/", type=str, help="Log folder path",required=True)
    parser.add_argument("--op_folder_path", default="output/", type=str, help="Output folde path",required=True)
    parser.add_argument("--ip_folder_path", default="input/", type=str, help="Input folder path",required=True)
    parser.add_argument("--root_run_dir", default=".", type=str, help="Root run directory",required=True)
    parser.add_argument("--resolution", default="3.6", type=str, help="Resolution")
    parser.add_argument("--job_id", default="XJFJKALJKLFJ3U471U", type=str, help="Job id")
    parser.add_argument("--input_map", default="input.mrc", type=str, help="Input map",required=True)
    parser.add_argument("--pdb_file_path", default="input.pdb", type=str, help="PDB file path",required=True)
    parser.add_argument("--pdb_name", default="UNKOWN", type=str, help="PDB name")
    parser.add_argument("--fasta_file_path", default="input.fasta", type=str, help="Fasta file path",required=True)
    parser.add_argument("--align_strategy", default="Manual alignment", type=str, help="Alignment strategy")
    parser.add_argument("--rosetta_path", default="rosetta", type=str, help="Rosetta path")
    args = parser.parse_args()
    # params = vars(args)
    return args