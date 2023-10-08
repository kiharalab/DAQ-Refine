# Based on the provided code snippet, let's incorporate the logic into the respective methods in the ProteinModeling class.
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

        # envrionment variables
        from sys import version_info
        self.python_version = f"{version_info.major}.{version_info.minor}"
        self.maxit_path = "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src/bin/maxit"
        self.RCSBROOT = "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src"

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
        self.rank_num = 1 #@param ["1", "2", "3", "4", "5"] {type:"raw"}
        self.color = '' #@param ["chain", "lDDT", "rainbow"]
        self.show_sidechains = False #@param {type:"boolean"}
        self.show_mainchains = False #@param {type:"boolean"}

        self.tag = ''
        self.jobname_prefix = ''
        self.pdb_filename = ''
        self.pdb_file = ''
        self.model_name = ''

        log_file = os.join(self.output_path, 'log.txt')
        logging.basicConfig(filename=log_file, level=logging.DEBUG)




    
    def print_parameters(self):
        # Extracted logic for printing help
        """Prints all the parameters of the instance."""
        logging.debug("=================================Parameters of the ProteinModeling instance:=================================")
        for attr, value in self.__dict__.items():
            logging.debug(f"{attr}: {value}")

        
    
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
                self.TrimDAQ(self.pdb_input_path, 0.0, self.input_path+'/1tmp.pdb')
            except Exception as e:
                logging.debug(f"Error while trimming DAQ: {e}")
                return False

            try:
                subprocess.run(["/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src/bin/maxit", "-input", self.input_path+"/1tmp.pdb", "-output", self.output_path+"/1tmp.cif", "-o", "1"], check=True)
                self.daq_file = self.output_path+'/1tmp.cif'
            except subprocess.CalledProcessError as e:
                logging.debug(f"Maxit subprocess failed: {e}")
                return False

            return True
        return False
    
    def prepare_trimmed_template(self):
        # Extracted logic for preparing trimmed template
        # Input target sequence
        # query_sequence = 'MGEVTAEEVEKFLDSNVSFAKQYYNLRYRAKVISDLLGPREAAVDFSNYHALNSVEESEIIFDLLRDFQDNLQAEKCVFNVMKKLCFLLQADRMSLFMYRARNGIAELATRLFNVHKDAVLEECLVAPDSEIVFPLDMGVVGHVALSKKIVNVPNTEEDEHFCDFVDTLTEYQTKNILASPIMNGKDVVAIIMVVNKVDGPHFTENDEEILLKYLNFANLIMKVFHLSYLHNCETRRGQILLWSGSKVFEELTDIERQFHKALYTVRAFLNCDRYSVGLLDMTKQKEFFDVWPVLMGEAPPYAGPRTPDGREINFYKVIDYILHGKEDIKVIPNPPPDHWALVSGLPTYVAQNGLICNIMNAPSEDFFAFQKEPLDESGWMIKNVLSMPIVNKKEEIVGVATFYNRKDGKPFDEMDETLMESLTQFLGWSVLNPDTYELMNKLENRKDIFQDMVKYHVKCDNEEIQTILKTREVYGKEPWECEEEELAEILQGELPDADKYEINKFHFSDLPLTELELVKCGIQMYYELKVVDKFHIPQEALVRFMYSLSKGYRRITYHNWRHGFNVGQTMFSLLVTGKLKRYFTDLEALAMVTAAFCHDIDHRGTNNLYQMKSQNPLAKLHGSSILERHHLEFGKTLLRDESLNIFQNLNRRQHEHAIHMMDIAIIATDLALYFKKRTMFQKIVDQSKTYETQQEWTQYMMLDQTRKEIVMAMMMTACDLSAITKPWEVQSKVALLVAAEFWEQGDLERTVLQQNPIPMMDRNKADELPKLQVGFIDFVCTFVYKEFSRFHEEITPMLDGITNNRKEWKALADEYETKMKGLEEEKQKQQAANQAAAGSQHGGKQPGGGPASKSCCVQ'
        # Remove whitespaces
        self.query_sequence = "".join(self.query_sequence.split())

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
                logging.debug("Could not find daq_file to rename.")
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
                logging.debug('Cannot find DAQ-score output file!!')
                exit(1)  # Or handle this case differently
            
            # print(f'User uploaded MSA file at {cust_msa_file}')
            
            try:
                a3m = self.ReadA3M(self.cust_msa_path)
                daq, good = self.ReadDAQ(self.pdb_input_path, 0.0, 0.0)
            except FileNotFoundError:
                logging.debug("MSA file or DAQ-score output file not found.")
                return False
            
            new_a3m = self.trim_a3m(a3m, daq, good)
            
            filename = os.path.join(self.output_path, 'trimmed_msa.a3m')
            self.save_a3m(filename, new_a3m)
            self.daq_msa = filename
            return self.daq_msa
        return
    
    def msa_settings(self):
        # Extracted logic for MSA settings
        self.msa_mode = "mmseqs2_uniref_env"
        self.pair_mode = "unpaired_paired"
        set_working_directory(self.output_path)

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
                    # print(line, end='')

                os.rename(self.custom_msa, self.a3m_file)
                self.queries_path=self.a3m_file
                logging.info(f"moving {self.custom_msa} to {self.a3m_file}")
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
                logging.info(f"moving {self.custom_msa} to {self.a3m_file}")
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

        set_working_directory('/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine')
        if not os.path.isfile("COLABFOLD_READY"):
            print("installing colabfold...")
            os.system("pip install -q --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/kiharalab/ColabFold'")
            os.system("pip install --upgrade dm-haiku")
            os.system("ln -s /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python{PYTHON_VERSION}/dist-packages/colabfold colabfold")
            os.system("ln -s /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine/lib/python{PYTHON_VERSION}/dist-packages/alphafold alphafold")
            # patch for jax > 0.3.25
            os.system("sed -i 's/weights = jax.nn.softmax(logits)/logits=jnp.clip(logits,-1e8,1e8);weights=jax.nn.softmax(logits)/g' alphafold/model/modules.py")
            os.system("touch COLABFOLD_READY")

        if USE_AMBER or USE_TEMPLATES:
            if not os.path.isfile("CONDA_READY"):
                print("installing conda...")
                os.system("wget -qnc https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh")
                os.system("bash Mambaforge-Linux-x86_64.sh -bfp /usr/local")
                os.system("conda config --set auto_update_conda false")
                os.system("touch CONDA_READY")

        if USE_TEMPLATES and not os.path.isfile("HH_READY") and USE_AMBER and not os.path.isfile("AMBER_READY"):
            print("installing hhsuite and amber...")
            os.system(f"conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 openmm=7.7.0 python='{PYTHON_VERSION}' pdbfixer")
            os.system("touch HH_READY")
            os.system("touch AMBER_READY")
        else:
            if USE_TEMPLATES and not os.path.isfile("HH_READY"):
                print("installing hhsuite...")
                os.system(f"conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 python='{PYTHON_VERSION}'")
                os.system("touch HH_READY")
            if USE_AMBER and not os.path.isfile("AMBER_READY"):
                print("installing amber...")
                os.system(f"conda install -y -c conda-forge openmm=7.7.0 python='{PYTHON_VERSION}' pdbfixer")
                os.system("touch AMBER_READY")

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
        
        check_gpu_with_torch()
        
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



        self.result_dir = f"{self.output_path}/results"
        os.makedirs(self.result_dir, exist_ok=True)
        self.log_filename = os.path.join(self.result_dir,"log.txt")
        if not os.path.isfile(self.log_filename) or 'logging_setup' not in globals():
            setup_logging(Path(self.log_filename))
        self.logging_setup = True

        self.queries, self.is_complex = get_queries(self.queries_path)
        self.model_type = set_model_type(self.is_complex, self.model_type)

        if "multimer" in self.model_type and self.max_msa is not None:
            self.use_cluster_profile = False
        else:
            self.use_cluster_profile = True

        logging.debug('=====================================DEBUG: PRINT PARAMETERS=====================================')
        self.print_parameters(self)
        logging.debug('=====================================DEBUG: PRINT PARAMETERS=====================================')

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
        results_zip = f"{self.jobname}.result.zip"
        os.system(f"zip -r {results_zip} {self.result_dir}")
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
    
    # TODO change the html show way to save figures
    def show_pdb(self,rank_num=1, show_sidechains=False, show_mainchains=False, color="lDDT"):
        self.model_name = f"rank_{rank_num}"
        view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
        view.addModel(open(self.pdb_file[0],'r').read(),'pdb')

        if color == "lDDT":
            view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
        elif color == "rainbow":
            view.setStyle({'cartoon': {'color':'spectrum'}})
        elif color == "chain":
            chains = len(self.queries[0][1]) + 1 if self.is_complex else 1
            for n,chain,color in zip(range(chains),alphabet_list,pymol_color_list):
                view.setStyle({'chain':chain},{'cartoon': {'color':color}})

        if show_sidechains:
            BB = ['C','O','N']
            view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
                                {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
            view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                                {'sphere':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
            view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                                {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
        if show_mainchains:
            BB = ['C','O','N','CA']
            view.addStyle({'atom':BB},{'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})

        view.zoomTo()
        return view
    
    def display_structure(self,results):
        self.rank_num = 1 #@param ["1", "2", "3", "4", "5"] {type:"raw"}
        self.color = "lDDT" #@param ["chain", "lDDT", "rainbow"]
        self.show_sidechains = False #@param {type:"boolean"}
        self.show_mainchains = False #@param {type:"boolean"}

        self.tag = results["rank"][0][self.rank_num - 1]
        self.jobname_prefix = ".custom" if self.msa_mode == "custom" else ""
        self.pdb_filename = f"{self.jobname}/{self.jobname}{self.jobname_prefix}_unrelaxed_{self.tag}.pdb"
        self.pdb_file = glob.glob(self.pdb_filename)

        

        self.show_pdb(self.rank_num, self.show_sidechains, self.show_mainchains, self.color).show()
        if self.color == "lDDT":
            plot_plddt_legend().show()

    def image_to_data_url(self,filename):
        ext = filename.split('.')[-1]
        prefix = f'data:image/{ext};base64,'
        with open(filename, 'rb') as f:
            img = f.read()
        return prefix + base64.b64encode(img).decode('utf-8')
    
    def display_html(self):

        # see: https://stackoverflow.com/a/53688522


        # pae = self.image_to_data_url(os.path.join(self.jobname,f"{self.jobname}{self.jobname_prefix}_pae.png"))
        # cov = self.image_to_data_url(os.path.join(self.jobname,f"{self.jobname}{self.jobname_prefix}_coverage.png"))
        # plddt = self.image_to_data_url(os.path.join(self.jobname,f"{self.jobname}{self.jobname_prefix}_plddt.png"))

        pae = self.save_image(os.path.join(self.jobname,f"{self.jobname}{self.jobname_prefix}_pae.png"))
        cov = self.save_image(os.path.join(self.jobname,f"{self.jobname}{self.jobname_prefix}_coverage.png"))
        plddt = self.save_image(os.path.join(self.jobname,f"{self.jobname}{self.jobname_prefix}_plddt.png"))
        display(HTML(f"""
        <style>
        img {{
            float:left;
        }}
        .full {{
            max-width:100%;
        }}
        .half {{
            max-width:50%;
        }}
        @media (max-width:640px) {{
            .half {{
            max-width:100%;
            }}
        }}
        </style>
        <div style="max-width:90%; padding:2em;">
        <h1>Plots for {escape(self.jobname)}</h1>
        <img src="{pae}" class="full" />
        <img src="{cov}" class="half" />
        <img src="{plddt}" class="half" />
        </div>
        """))
    
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

        print("INFO: STEP-1 Input Protein Sequence and DAQ result file started")
        # logging.info("INFO: STEP-1 Input Protein Sequence and DAQ result file started")

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

        print("INFO: STEP-1 Input Protein Sequence and DAQ result file Done")
        print("INFO: STEP-2 Modeling Part started")

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

        try:
            results = self.prediction()
        except Exception as e:
            print(f"Error in prediction(): {e}")
            exit(1)

        print("Prediction finished.")
        print("INFO: STEP-2 Modeling Part Done")
        print("INFO: STEP-3 Save Results Started")
        # self.dispaly_structure(results)
        print("INFO: STEP-3 Save Results Done")


        print("======================DAQ-Refine finished==================")
        print("Selected strategy mode:", self.str_mode)
        print("Input query sequence:", self.query_sequence)
        print("Job ID:", self.jobname)
        print("Number of models to use:", self.num_relax)
        print("Template mode:", self.template_mode)
        print("Output directory:", self.output_path)



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

def get_arguments():
    parser = argparse.ArgumentParser(description='STEP-1: Input Protein Sequence and DAQ result file')

    # Add arguments
    parser.add_argument('--str_mode', type=str, default='strategy 2',
                        help='Select the DAQ-refine strategy. Choices are Vanilla AF2, strategy 1, and strategy 2.',required=True)
    
    parser.add_argument('--query_sequence', type=str, default='MGEVTAEEVEKFLDSNVSFAKQYYNLRYRAKVISDLLGPREAAVDFSNYHALNSVEESEIIFDLLRDFQDNLQAEKCVFNVMKKLCFLLQADRMSLFMYRARNGIAELATRLFNVHKDAVLEECLVAPDSEIVFPLDMGVVGHVALSKKIVNVPNTEEDEHFCDFVDTLTEYQTKNILASPIMNGKDVVAIIMVVNKVDGPHFTENDEEILLKYLNFANLIMKVFHLSYLHNCETRRGQILLWSGSKVFEELTDIERQFHKALYTVRAFLNCDRYSVGLLDMTKQKEFFDVWPVLMGEAPPYAGPRTPDGREINFYKVIDYILHGKEDIKVIPNPPPDHWALVSGLPTYVAQNGLICNIMNAPSEDFFAFQKEPLDESGWMIKNVLSMPIVNKKEEIVGVATFYNRKDGKPFDEMDETLMESLTQFLGWSVLNPDTYELMNKLENRKDIFQDMVKYHVKCDNEEIQTILKTREVYGKEPWECEEEELAEILQGELPDADKYEINKFHFSDLPLTELELVKCGIQMYYELKVVDKFHIPQEALVRFMYSLSKGYRRITYHNWRHGFNVGQTMFSLLVTGKLKRYFTDLEALAMVTAAFCHDIDHRGTNNLYQMKSQNPLAKLHGSSILERHHLEFGKTLLRDESLNIFQNLNRRQHEHAIHMMDIAIIATDLALYFKKRTMFQKIVDQSKTYETQQEWTQYMMLDQTRKEIVMAMMMTACDLSAITKPWEVQSKVALLVAAEFWEQGDLERTVLQQNPIPMMDRNKADELPKLQVGFIDFVCTFVYKEFSRFHEEITPMLDGITNNRKEWKALADEYETKMKGLEEEKQKQQAANQAAAGSQHGGKQPGGGPASKSCCVQ',
                        help='Input target protein sequence.')

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

    args = parser.parse_args()

    return args

def main():
    # Get arguments (this function needs to be implemented based on the original code)
    args = get_arguments()
    
    # Create an instance of the ProteinModeling class
    modeling = Daqrefine(
        args=args
    )
    
    # Run the modeling process
    modeling.run_modeling()

if __name__ == '__main__':
    main()

