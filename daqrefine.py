
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

def check_gpu_with_torch():
    # Check if GPU is available
    if torch.cuda.is_available():
        print("GPU is available!")
        # Get the name of the GPU
        print(f"GPU Name: {torch.cuda.get_device_name(0)}")
        # Get memory information about the GPU
        total_memory = torch.cuda.get_device_properties(0).total_memory
        allocated_memory = torch.cuda.memory_allocated()
        cached_memory = torch.cuda.memory_cached()
        free_memory = total_memory - allocated_memory - cached_memory
        print(f"Total Memory: {total_memory / (1024 ** 2)} MB")
        print(f"Allocated Memory: {allocated_memory / (1024 ** 2)} MB")
        print(f"Cached Memory: {cached_memory / (1024 ** 2)} MB")
        print(f"Free Memory: {free_memory / (1024 ** 2)} MB")
    else:
        print("GPU is not available.")

def set_working_directory(path):
    """
    Set the current working directory to the specified path.
    
    Args:
    - path (str): Path to set as the current working directory.
    """
    try:
        os.chdir(path)
        print(f"Current working directory set to {path}")
    except Exception as e:
        print(f"Error setting working directory: {e}")

class Daqrefine:
    def __init__(self,args):


        # Initialize class variables

        # pass args to the class
        self.cust_msa_path = args.cust_msa_path
        self.output_path = args.output_path
        self.input_path = args.input_path
        self.query_sequence = args.query_sequence
        self.jobname = args.jobname        
        self.num_relax = args.num_relax
        self.template_mode = args.template_mode
        self.str_mode = args.str_mode
        self.pdb_input_path = args.pdb_input_path
        self.resolution = args.resolution

        # envrionment variables
        from sys import version_info
        self.python_version = f"{version_info.major}.{version_info.minor}"
        self.emweb_path = "/bio/kihara-web/www/em/emweb-jobscheduler"
        self.emweb_daqrefine_path = os.path.join(self.emweb_path,"algorithms/DAQ-Refine")
        self.emweb_daq_path = os.path.join(self.emweb_path,"algorithms/DAQ")
        self.mmalign_path = os.path.join(self.emweb_daqrefine_path,"MMalign")
        
        self.RCSBROOT = os.path.join(self.emweb_daqrefine_path,"maxit-v11.100-prod-src")
        self.maxit_path = os.path.join(self.RCSBROOT,"bin/maxit")

        # initialize these parameters in step-1
        self.daq_file = ''
        self.daq_msa = ''
        self.use_templates = False
        self.use_amber = False
        self.dispaly_images = True
        self.queries_path = ''

        self.msa_mode = ''
        self.pair_mode = ''
        self.a3m_file = ''
        self.custom_msa = ''

        # parameters in modeling, initialize these parameters in step-2
        self.model_type = ""
        self.num_recycles = "0"
        self.recycle_early_stop_tolerance = ""
        self.pairing_strategy = ""

        self.max_msa = "" #@param ["auto", "512:1024", "256:512", "64:128", "32:64", "16:32"]
        self.num_seeds = 1 #@param [1,2,4,8,16] {type:"raw"}
        self.use_dropout = False #@param {type:"boolean"}

        self.num_recycles = ''
        self.recycle_early_stop_tolerance = None 
        self.save_all = False #@param {type:"boolean"}
        self.save_recycles = False #@param {type:"boolean"}
        self.save_to_google_drive = False #@param {type:"boolean"}
        self.dpi = 200 #@param {type:"integer"}

        self.K80_chk = ''
        self.logging_setup = False
        self.queries = ''
        self.is_complex = False
        self.result_dir = ''

        # parameters used in display
        self.dispaly_images = False
        self.color = '' #@param ["chain", "lDDT", "rainbow"]
        self.show_sidechains = False #@param {type:"boolean"}
        self.show_mainchains = False #@param {type:"boolean"}

        self.tag = ''
        self.jobname_prefix = ''
        self.pdb_filename = ''
        # self.pdb_file = ''
        self.model_name = ''
        

        # parrameters in rerun DAQ
        self.rerun_daq_result_path = ''
        self.final_pdb_path = ''

        log_file = os.path.join(self.output_path, 'log.txt')
        # logging.basicConfig(filename=log_file, level=DEBUG)




    def clean_up(self,args):
        set_working_directory(args.output_path)
        
        self._cleanup_folders(["S1_results","S2_results", "Vanilla_AF2_results", "DAQ", "template"],args)
        self._cleanup_files_with_suffixes([".a3m", ".csv", ".txt"],args)

    def _cleanup_folders(self, folders,args):
        for folder in folders:
            folder_path = os.path.join(args.output_path, folder)
            if os.path.exists(folder_path) and os.path.isdir(folder_path):
                try:
                    shutil.rmtree(folder_path)
                    print(f"'{folder}' folder has been removed.")
                except Exception as e:
                    print(f"Error deleting '{folder}' folder. Reason: {e}")
            else:
                print(f"'{folder}' folder does not exist in the current directory.")

    def _cleanup_files_with_suffixes(self, suffixes,args):
        for filename in os.listdir(args.output_path):
            if any(filename.endswith(end) for end in suffixes):
                file_path = os.path.join(args.output_path, filename)
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Error deleting {file_path}. Reason: {e}")


    
    def print_parameters(self):
        # Extracted logic for printing help
        """Prints all the parameters of the instance."""
        print("=================================Parameters of the DaqRefine instance:=================================")
        for attr, value in self.__dict__.items():
            print(f"{attr}: {value}")



    def pdb_to_fasta(self,pdb_path):
        # Create a PDB parser object
        parser = PDBParser()

        # Parse the PDB file
        structure = parser.get_structure("pdb_structure", pdb_path)

        # Create a sequence object to store the protein sequence
        seq = ""

        # Loop over all the chains in the protein
        for model in structure:
            for chain in model:
                # Create a string to store the amino acid sequence for this chain
                chain_seq = ""

                # Loop over all the residues in the chain
                for residue in chain:
                    # Check if the residue is an amino acid
                    if residue.get_resname() not in ["HOH", "HOH2", "UNX", "UNL", "UNX1"]:
                        # Get the one-letter amino acid code for the residue
                        aa = residue.get_resname().upper()

                        # Add the amino acid code to the chain sequence string
                        chain_seq += aa

                # Create a SeqRecord object for this chain sequence
                record = SeqRecord(Seq(chain_seq), id=chain.get_id(), description="")

                # Add the sequence to the overall protein sequence
                seq += str(record.seq)

        # Create a SeqRecord object for the protein sequence
        seq = seq.replace(" ","")
        protein_record = SeqRecord(Seq(seq), id=structure.id, description="")

        return seq
        # Return the protein sequence in FASTA format
        # return protein_record.format("fasta")

# import sys
# fasta_string = pdb_to_fasta(sys.argv[1])
# print(fasta_string)  
    
    def TrimDAQ(self, filename, cutoff, outfile):
        daq = []
        PDB = {}
        lines = ''
        with open(filename) as f:
            for li in f:
                if li.startswith('ATOM'):
                    li = li.strip()
                    resn = int(li[23:27])
                    sco = float(li[60:66])
                    x = float(li[30:38])
                    y = float(li[38:46])
                    z = float(li[46:54])
                    PDB[str(resn)] = [x, y, z, sco]
                    if sco < cutoff:
                        daq.append(resn)
                    else:
                        lines = lines + li + '\n'
        
        with open(outfile, 'w') as out:
            out.write(lines)
    
    def add_hash(self,x, y):
        return x + "_" + hashlib.sha1(y.encode()).hexdigest()[:5]


    def get_input(self):
        # Extracted logic for getting inputs

        if self.str_mode in ("strategy 1", "strategy 2"):
            try:
                self.TrimDAQ(self.pdb_input_path, 0.0, self.output_path+'/1tmp.pdb')
            except Exception as e:
                print(f"Error while trimming DAQ: {e}")
                return False

            try:
                subprocess.run(["/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src/bin/maxit", "-input", self.output_path+"/1tmp.pdb", "-output", self.output_path+"/1tmp.cif", "-o", "1"], check=True)
                self.daq_file = self.output_path+'/1tmp.cif'
            except subprocess.CalledProcessError as e:
                print(f"Maxit subprocess failed: {e}")
                return False

            return True
        return True
    
    def prepare_trimmed_template(self):
        # remove all whitespaces
        self.query_sequence = "".join(self.query_sequence.split())
        # remove all \r, \n
        self.query_sequence = self.query_sequence.replace('\r', '').replace('\n', '')

        # Job name
        self.basejobname = "".join(self.jobname.split())
        self.basejobname = re.sub(r'\W+', '', self.basejobname)
        self.jobname = self.add_hash(self.basejobname, self.query_sequence)

        set_working_directory(self.output_path)

        while os.path.isfile(f"{self.jobname}.csv"):
            self.jobname = self.add_hash(self.basejobname, ''.join(random.sample(self.query_sequence, len(self.query_sequence))))

        with open(f"{self.jobname}.csv", "w") as text_file:
            text_file.write(f"id,sequence\n{self.jobname},{self.query_sequence}")

        self.queries_path = f"{self.jobname}.csv"

        # Number of models to use
        self.num_relax = int(self.num_relax)
        self.use_amber = self.num_relax > 0

        if self.str_mode == "strategy 1" or self.str_mode == "strategy 2":
            self.custom_template_path = f"{self.output_path}/template"
            os.makedirs(self.custom_template_path, exist_ok=True)
                        
            self.use_templates = True
                
            # Assume daq_file is provided as a local path by the user
            
            try:
                os.rename(self.daq_file, f"template/1tmp.cif")
            except FileNotFoundError:
                print("Could not find daq_file to rename.")
                return False

            #     return True
                
            self.template_mode = "custom"

        elif self.template_mode == "pdb70":
            self.use_templates = True
            self.custom_template_path = None

        elif self.template_mode == "custom":
            self.custom_template_path = f"{self.output_path}/template"
            # TODO add logic later
            
            self.use_templates = True

        else:
            self.custom_template_path = None
            self.use_templates = False
        return
    
    def trim_a3m(self,a3m,daq,good):

        out=[]
        for ith in range(len(a3m)):
            name,seq = a3m[ith]
            new_seq=''
            if ith == 0:#query
                new_seq = seq
            else:
                pos=0
                for aa in seq:
                    if aa == aa.lower() and aa != '-':
                        continue
                    pos = pos + 1
                    if pos in daq: #selected bad regions or missing regions
                        new_seq = new_seq + aa
                    elif not pos in good:
                        new_seq = new_seq + aa
                    elif aa == aa.upper():
                        new_seq = new_seq + '-'
            out.append([name,new_seq])
        return out

    def ReadA3M(self,filename):
        A3M=[]
        with open(filename) as f:
            for li in f:
                if li.startswith('#'):
                    continue
                if li.startswith('>'):
                    name = li.strip()
                else:
                    seq = li.strip()
                    A3M.append([name,seq])
        return A3M
    
    def ReadDAQ(self,filename,cutoff,dist_cut):
        daq=[]
        PDB={}
        with open(filename) as f:
            for li in f:
                #print(li)
                if li.startswith('ATOM') and li[13:16]=='CA ':
                    li=li.strip()
                    resn = int(li[23:27])
                    sco = float(li[60:66])
                    x = float(li[30:38])
                    y = float(li[38:46])
                    z = float(li[46:54])
                    PDB[str(resn)]=[x,y,z,sco]
                    if sco < cutoff:
                        #print(sco)
                        daq.append(resn)
        # print(daq)
        daq2=[]
        for resn in PDB:
            if int(resn) in daq:
                continue
            #Distance check
            x1=PDB[resn][0]
            y1=PDB[resn][1]
            z1=PDB[resn][2]
            for resn2 in PDB:
                if not int(resn2) in daq:
                    continue
                x2=PDB[resn2][0]
                y2=PDB[resn2][1]
                z2=PDB[resn2][2]
                dist=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                if dist <=dist_cut*dist_cut:#close
                    daq2.append(int(resn))
                    break
        #print('LowDAQ',daq,'Extended',daq2)
        daq=list(set(daq+daq2))
        #Others
        goodpos=[]
        for resn in PDB:
            if int(resn) in daq:
                continue
            goodpos.append(int(resn))
        # print('HighDAQ',goodpos)

        return daq,goodpos
    
    def save_a3m(self,file,a3m):
        lines=''
        for name,seq in a3m:
            lines = lines + name+'\n'
            lines = lines + seq+'\n'
        with open(file,'w') as out:
            out.write(lines)

    def msa(self):
        # Extracted logic for MSA
        # pdb_input_path = args.pdb_input_path
        if self.str_mode == "strategy 2":

            if not os.path.isfile(self.pdb_input_path):
                print('Cannot find DAQ-score output file!!')
                exit(1)  # Or handle this case differently
            
            # print(f'User uploaded MSA file at {cust_msa_file}')
            
            # print(self.cust_msa_path)
            try:
                a3m = self.ReadA3M(self.cust_msa_path)
                daq, good = self.ReadDAQ(self.pdb_input_path, 0.0, 0.0)
            except FileNotFoundError:
                print("MSA file not found.")
                exit(1)
            
            new_a3m = self.trim_a3m(a3m, daq, good)
            
            filename = os.path.join(self.output_path, 'trimmed_msa.a3m')
            self.save_a3m(filename, new_a3m)
            self.daq_msa = filename
            # file_path = self.daq_msa
            # file_size = os.path.getsize(file_path)

            # print(f"The size of '{file_path}' is {file_size} bytes.")
            return self.daq_msa
        return
    
    def msa_settings(self):
        # Extracted logic for MSA settings
        self.msa_mode = "mmseqs2_uniref_env"
        self.pair_mode = "unpaired_paired"
        set_working_directory(self.output_path)
        # file_path = self.daq_msa
        # file_size = os.path.getsize(file_path)

        # print(f"The size of '{file_path}' is {file_size} bytes.")

        if self.str_mode == "strategy 2":
            self.msa_mode = "custom"
            self.a3m_file = f"{self.jobname}.custom.a3m"
            if not os.path.isfile(self.a3m_file):
                self.custom_msa = self.daq_msa
                header = 0
                # import fileinput
                for line in fileinput.FileInput(self.custom_msa,inplace=1):
                    if line.startswith(">"):
                        header = header + 1
                    if not line.rstrip():
                        continue
                    if line.startswith(">") == False and header == 1:
                        self.query_sequence = line.rstrip()
                    print(line, end='')

                os.rename(self.custom_msa, self.a3m_file)
                self.queries_path=self.a3m_file
                print(f"moving {self.custom_msa} to {self.a3m_file}")
                file_path = self.a3m_file
                file_size = os.path.getsize(file_path)

                print(f"The size of '{file_path}' is {file_size} bytes.")
        elif self.msa_mode.startswith("mmseqs2"):
            self.a3m_file = f"{self.jobname}.a3m"
        elif self.msa_mode == "custom":
            self.a3m_file = f"{self.jobname}.custom.a3m"
            if not os.path.isfile(self.a3m_file):
                # self.custom_msa = input("Enter path to the custom MSA file: ")
                # TODO: remove this logic in the future
                header = 0
                for line in fileinput.FileInput(self.custom_msa, inplace=1):
                    if line.startswith(">"):
                        header += 1
                    if not line.rstrip():
                        continue
                    if not line.startswith(">") and header == 1:
                        query_sequence = line.rstrip()
                    # print(line, end='')
                os.rename(self.custom_msa, self.a3m_file)
                self.queries_path = self.a3m_file
                print(f"moving {self.custom_msa} to {self.a3m_file}")
        else:
            self.a3m_file = f"{self.jobname}.single_sequence.a3m"
            with open(self.a3m_file, "w") as text_file:
                # query_sequence = input("Enter query sequence: ")
                text_file.write(">1\n%s" % self.query_sequence)
    
    def advanced_setting(self):
        # Extracted logic for advanced settings
        self.model_type = "auto"
        self.num_recycles = "1"
        self.recycle_early_stop_tolerance = "auto"
        self.pairing_strategy = "greedy"

        self.max_msa = "auto" #@param ["auto", "512:1024", "256:512", "64:128", "32:64", "16:32"]
        self.num_seeds = 1 #@param [1,2,4,8,16] {type:"raw"}
        self.use_dropout = False #@param {type:"boolean"}

        self.num_recycles = None if self.num_recycles == "auto" else int(self.num_recycles)
        self.recycle_early_stop_tolerance = None if self.recycle_early_stop_tolerance == "auto" else float(self.recycle_early_stop_tolerance)
        if self.max_msa == "auto": self.max_msa = None

        #@markdown #### Save settings
        self.save_all = False #@param {type:"boolean"}
        self.save_recycles = False #@param {type:"boolean"}
        self.save_to_google_drive = False #@param {type:"boolean"}
        #@markdown -  if the save_to_google_drive option was selected, the result zip will be uploaded to your Google Drive
        self.dpi = 200 #@param {type:"integer"}
        #@markdown - set dpi for image resolution

        # if save_to_google_drive:
        #     from pydrive.drive import GoogleDrive
        #     from pydrive.auth import GoogleAuth
        #     from google.colab import auth
        #     from oauth2client.client import GoogleCredentials
        #     auth.authenticate_user()
        #     gauth = GoogleAuth()
        #     gauth.credentials = GoogleCredentials.get_application_default()
        #     drive = GoogleDrive(gauth)
        #     print("You are logged into Google Drive and are good to go!")
    
    def check_dependencies(self):
        USE_AMBER = self.use_amber
        USE_TEMPLATES = self.use_templates
        PYTHON_VERSION = self.python_version

        def is_python_module_installed(module_name):
            try:
                subprocess.check_call([f"python{PYTHON_VERSION}", "-c", f"import {module_name}"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                return True
            except subprocess.CalledProcessError:
                return False
        
        def is_conda_package_installed(package_name):
            try:
                result = subprocess.check_output(["conda", "list", package_name], stderr=subprocess.STDOUT, text=True)
                if package_name in result:
                    return True
                else:
                    return False
            except subprocess.CalledProcessError:
                return False

        def is_conda_installed():
            try:
                subprocess.check_call(["conda", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                return True
            except subprocess.CalledProcessError:
                return False

        set_working_directory('/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine')

        if not is_python_module_installed("colabfold"):
            print("installing colabfold...")
            os.system("pip install -q --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/kiharalab/ColabFold'")
            os.system("pip install --upgrade dm-haiku")
            os.system(f"ln -s /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python{PYTHON_VERSION}/dist-packages/colabfold colabfold")
            os.system(f"ln -s /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python{PYTHON_VERSION}/dist-packages/alphafold alphafold")
            # patch for jax > 0.3.25
            os.system("sed -i 's/weights = jax.nn.softmax(logits)/logits=jnp.clip(logits,-1e8,1e8);weights=jax.nn.softmax(logits)/g' alphafold/model/modules.py")

        if not is_conda_installed():
            print("installing conda...")
            os.system("wget -qnc https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh")
            os.system("bash Miniconda3-latest-Linux-x86_64.sh -bfp /usr/local")
            os.system("conda config --set auto_update_conda false")

        if USE_TEMPLATES and not is_conda_package_installed("hhsuite") and USE_AMBER and not is_conda_package_installed("openmm"):
            print("installing hhsuite and amber...")
            os.system(f"conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 openmm=7.7.0 python='{PYTHON_VERSION}' pdbfixer")
        else:
            if USE_TEMPLATES and not is_conda_package_installed("hhsuite"):
                print("installing hhsuite...")
                os.system(f"conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 python='{PYTHON_VERSION}'")
            if USE_AMBER and not is_conda_package_installed("openmm"):
                print("installing amber...")
                os.system(f"conda install -y -c conda-forge openmm=7.7.0 python='{PYTHON_VERSION}' pdbfixer")

    def input_features_callback(self,input_features):
        if self.display_images:
            plot_msa_v2(input_features)
            plt.show()
            plt.close()

    def prediction_callback(self,protein_obj, length,prediction_result, input_features, mode):
        self.model_name, relaxed = mode

        if not relaxed:
            if self.display_images:
                fig = plot_protein(protein_obj, Ls=length, dpi=150)
                plt.show()
                plt.close()
    
    def prediction(self):
        # Extracted logic for prediction
        warnings.simplefilter(action='ignore', category=FutureWarning)
        warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)
        self.display_images = False #@param {type:"boolean"}
        set_working_directory(self.output_path)
        
        
        
        try:
            self.K80_chk = os.popen('nvidia-smi | grep "Tesla K80" | wc -l').read()
        except:
            self.K80_chk = "0"
        pass
        if "1" in self.K80_chk:
            print("WARNING: found GPU Tesla K80: limited to total length < 1000")
        if "TF_FORCE_UNIFIED_MEMORY" in os.environ:
            del os.environ["TF_FORCE_UNIFIED_MEMORY"]
        if "XLA_PYTHON_CLIENT_MEM_FRACTION" in os.environ:
            del os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"]


        # For some reason we need that to get pdbfixer to import
        # /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python3.8/site-packages
        if self.use_amber and f"/bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python{self.python_version}/site-packages/" not in sys.path:
            sys.path.insert(0, f"/bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python{self.python_version}/site-packages/")



        if self.str_mode == "Vanilla AF2":
            self.result_dir = f"{self.output_path}/Vanilla_AF2_results"
        elif self.str_mode ==  "strategy 1":
            self.result_dir = f"{self.output_path}/S1_results"
        else:
            self.result_dir = f"{self.output_path}/S2_results"
        os.makedirs(self.result_dir, exist_ok=True)
        self.log_filename = os.path.join(self.result_dir,"colabfold_log.txt")
        if not os.path.isfile(self.log_filename) or 'logging_setup' not in globals():
            setup_logging(Path(self.log_filename))
        self.logging_setup = True

        self.queries, self.is_complex = get_queries(self.queries_path)
        self.model_type = set_model_type(self.is_complex, self.model_type)

        if "multimer" in self.model_type and self.max_msa is not None:
            self.use_cluster_profile = False
        else:
            self.use_cluster_profile = True

        # print('=====================================DEBUG: PRINT PARAMETERS=====================================')
        # self.print_parameters()
        # print('=====================================DEBUG: PRINT PARAMETERS=====================================')

        download_alphafold_params(self.model_type, Path("."))
        results = run(
            queries=self.queries,
            result_dir=self.result_dir,
            use_templates=self.use_templates,
            custom_template_path=self.custom_template_path,
            num_relax=self.num_relax,
            msa_mode=self.msa_mode,
            model_type=self.model_type,
            num_models=5,
            num_recycles=self.num_recycles,
            recycle_early_stop_tolerance=self.recycle_early_stop_tolerance,
            num_seeds=self.num_seeds,
            use_dropout=self.use_dropout,
            model_order=[1,2,3,4,5],
            is_complex=self.is_complex,
            data_dir=Path("."),
            keep_existing_results=False,
            rank_by="auto",
            pair_mode=self.pair_mode,
            pairing_strategy=self.pairing_strategy,
            stop_at_score=float(100),
            prediction_callback=self.prediction_callback,
            dpi=self.dpi,
            zip_results=False,
            save_all=self.save_all,
            max_msa=self.max_msa,
            use_cluster_profile=self.use_cluster_profile,
            input_features_callback=self.input_features_callback,
            save_recycles=self.save_recycles,
        )
        # results_zip = f"{self.jobname}.result.zip"
        # os.system(f"zip -r {results_zip} {self.result_dir}")
        return results
    
    def save_image(self, filename):
        # Define the storage location
        storage_location = f"{self.output_path}\images"
        
        # Create the directory if it doesn't exist
        if not os.path.exists(storage_location):
            os.makedirs(storage_location)
        
        # Copy the image to the storage location
        destination_path = os.path.join(storage_location, os.path.basename(filename))
        shutil.copyfile(filename, destination_path)
        
        return destination_path
    
    
    
    def align_structure(self,results):
        self.color = "lDDT" #@param ["chain", "lDDT", "rainbow"]
        self.show_sidechains = False #@param {type:"boolean"}
        self.show_mainchains = False #@param {type:"boolean"}

            
        self.rerun_daq_result_path = os.path.join(self.output_path,"DAQ")
        # input_map = os.path.join(self.output_path,"input_resize.mrc")
        # input_pdb = os.path.join(self.rerun_daq_result_path,"DAQ/input.pdb")

        set_working_directory(self.output_path)

        # mkdir DAQ
        if not os.path.exists(self.rerun_daq_result_path):
            os.mkdir(self.rerun_daq_result_path)

        rank_range = []
        if self.str_mode == "strategy 1":
            rank_range = [1,2]
        else:
            rank_range = [1,2,3,4,5]
        for rank_num in rank_range:
            self.tag = results["rank"][0][rank_num - 1]
            self.jobname_prefix = ".custom" if self.msa_mode == "custom" else ""
            self.pdb_filename = f"{self.result_dir}/{self.jobname}{self.jobname_prefix}_unrelaxed_{self.tag}.pdb"



            if not os.path.exists(self.pdb_filename):
                print(f"File '{self.pdb_filename}' not found.")
                exit(0)

            if self.str_mode=="strategy 1":
                temp_tag = self.tag + "_s1"
            elif self.str_mode=="strategy 2":
                temp_tag = self.tag + "_s2"
            else:
                temp_tag = self.tag + "_af2"
            # MMalign
            self.mmalign(self.pdb_filename,temp_tag)
            # ROSETTA3 relaxation
            self.rosetta_relaxation(os.path.join(self.output_path,f"DAQ/{temp_tag}"))

    
    def mmalign(self,input_pdb,tag):
        set_working_directory(self.output_path)
        if not os.path.exists(os.path.join(self.output_path,f"DAQ/{tag}")):
            os.mkdir(os.path.join(self.output_path,f"DAQ/{tag}"))
        aligned_output_dir = os.path.join(self.output_path,f"DAQ/{tag}")
        try:
            print('align structure to input map')
            result = subprocess.run([self.mmalign_path, input_pdb, self.pdb_input_path, "-o", f"{aligned_output_dir}/input.pdb"], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing {self.mmalign_path}. Return code: {e.returncode}")
            print(f"Output: {e.output}")
            print(f"Error: {e.stderr}")
        except Exception as e:
            print(f"Unexpected error: {e}")
    
    # TODO figure out the parameters for paths
    def rosetta_relaxation(self,input_dir):
        set_working_directory(input_dir)
        shutil.copy(os.path.join(self.emweb_daqrefine_path, "B_relax_density.xml"),"B_relax_density.xml")
        shutil.copy(os.path.join(self.output_path,"input_resize.mrc"),"MAP.mrc")
        ROSETTA3 = os.path.expanduser("/home/kihara/gterashi/bin/rosetta_bin_linux_2018.33.60351_bundle/main")
        command = [
            os.path.join(ROSETTA3, "source/bin/rosetta_scripts.static.linuxgccrelease"),
            "-database", os.path.join(ROSETTA3, "database/"),
            "-in::file::s", "input.pdb",
            "-parser::protocol", "B_relax_density.xml",
            "-ignore_unrecognized_res",
            "-edensity::mapreso", self.resolution,
            "-edensity::cryoem_scatterers",
            "-crystal_refine",
            "-beta",
            "-out::suffix", "_relax",
            "-default_max_cycles", "200"
        ]

        result = subprocess.run(command, capture_output=True, text=True)
        
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        return       


    def run_modeling(self):
        # Extracted logic for running the modeling process
        try:
            # Set environment variables
            os.environ['RCSBROOT'] = self.RCSBROOT
            os.environ['PATH'] += ":RCSBROOT$/bin"
        except Exception as e:
            print(f"Error setting environment variables: {e}")
            exit(1)

        try:
            input_success = self.get_input()
            if not input_success:
                print("Exiting due to error in input.")
                exit(1)
        except Exception as e:
            print(f"Error in get_input(): {e}")
            exit(1)

        print("INFO: Input Protein Sequence and DAQ result file started")
        # print("INFO: STEP-1 Input Protein Sequence and DAQ result file started")

        try:
            self.prepare_trimmed_template()
            print("Prepare trimmed template finished.")
        except Exception as e:
            print(f"Error in prepare_trimmed_template(): {e}")
            exit(1)

        try:
            self.msa()
            print("MSA finished.(if applicable)")
        except Exception as e:
            print(f"Error in msa(): {e}")
            exit(1)

        print("INFO: Input Protein Sequence and DAQ result file Done")
        print("INFO: Modeling Part started")

        try:
            self.msa_settings()
            print("MSA settings finished.")
        except Exception as e:
            print(f"Error in msa_settings(): {e}")
            exit(1)

        try:
            self.advanced_setting()
            print("Advanced settings finished.")
        except Exception as e:
            print(f"Error in advanced_setting(): {e}")
            exit(1)

        try:
            self.check_dependencies()
            print("Install dependencies finished.")
        except Exception as e:
            print(f"Error in check_dependencies(): {e}")
            exit(1)

        results = None
        try:
            
            check_gpu_with_torch()
            results = self.prediction()
            print("Prediction finished.")
        except Exception as e:
            print(f"Error in prediction(): {e}")
            exit(1)


        print("INFO: Modeling Part Done")
        try:
            self.align_structure(results)
            print("Align structures finished.")
        except Exception as e:
            print(f"Error in align_structure(results): {e}")
            exit(1)
        


        print("====================================================================Modeling finished====================================================================")
        print("Selected strategy mode:", self.str_mode)
        print("Input query sequence:", self.query_sequence)
        print("Job ID:", self.jobname)
        print("Number of models to use:", self.num_relax)
        print("Template mode:", self.template_mode)
        print("Output directory:", self.rerun_daq_result_path)
        print("====================================================================Modeling finished====================================================================")
