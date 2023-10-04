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


from sys import version_info
python_version = f"{version_info.major}.{version_info.minor}"
maxit_path = "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src/bin/maxit"
RCSBROOT = "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src"
# STEP-1: Input Protein Sequence and DAQ result file

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

def TrimDAQ(filename,cutoff,outfile):
    daq=[]
    PDB={}
    lines=''
    with open(filename) as f:
        for li in f:
            #print(li)
            if li.startswith('ATOM'):
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
                else:
                    lines=lines+li+'\n'
    #print(lines)
    with open(outfile,'w') as out:
        out.write(lines)

daq_file = ''
daq_msa = ''
template_mode= ''
query_sequence = ''
use_templates = False
use_amber = False
dispaly_images = True
jobname=''
queries_path = ''
def add_hash(x, y):
    return x + "_" + hashlib.sha1(y.encode()).hexdigest()[:5]

def get_input(args):
    global daq_file

    if args.str_mode in ("strategy 1", "strategy 2"):
        try:
            TrimDAQ(args.pdb_input_path, 0.0, args.input_path+'/1tmp.pdb')
        except Exception as e:
            print(f"Error while trimming DAQ: {e}")
            return False

        try:
            subprocess.run(["/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine/maxit-v11.100-prod-src/bin/maxit", "-input", args.input_path+"/1tmp.pdb", "-output", args.output_path+"/1tmp.cif", "-o", "1"], check=True)
            daq_file = args.output_path+'/1tmp.cif'
        except subprocess.CalledProcessError as e:
            print(f"Maxit subprocess failed: {e}")
            return False

        return True
    return False



def prepare_trimmed_template(args):
    global daq_file, daq_msa, template_mode, query_sequence
    # Input target sequence
    query_sequence = 'MGEVTAEEVEKFLDSNVSFAKQYYNLRYRAKVISDLLGPREAAVDFSNYHALNSVEESEIIFDLLRDFQDNLQAEKCVFNVMKKLCFLLQADRMSLFMYRARNGIAELATRLFNVHKDAVLEECLVAPDSEIVFPLDMGVVGHVALSKKIVNVPNTEEDEHFCDFVDTLTEYQTKNILASPIMNGKDVVAIIMVVNKVDGPHFTENDEEILLKYLNFANLIMKVFHLSYLHNCETRRGQILLWSGSKVFEELTDIERQFHKALYTVRAFLNCDRYSVGLLDMTKQKEFFDVWPVLMGEAPPYAGPRTPDGREINFYKVIDYILHGKEDIKVIPNPPPDHWALVSGLPTYVAQNGLICNIMNAPSEDFFAFQKEPLDESGWMIKNVLSMPIVNKKEEIVGVATFYNRKDGKPFDEMDETLMESLTQFLGWSVLNPDTYELMNKLENRKDIFQDMVKYHVKCDNEEIQTILKTREVYGKEPWECEEEELAEILQGELPDADKYEINKFHFSDLPLTELELVKCGIQMYYELKVVDKFHIPQEALVRFMYSLSKGYRRITYHNWRHGFNVGQTMFSLLVTGKLKRYFTDLEALAMVTAAFCHDIDHRGTNNLYQMKSQNPLAKLHGSSILERHHLEFGKTLLRDESLNIFQNLNRRQHEHAIHMMDIAIIATDLALYFKKRTMFQKIVDQSKTYETQQEWTQYMMLDQTRKEIVMAMMMTACDLSAITKPWEVQSKVALLVAAEFWEQGDLERTVLQQNPIPMMDRNKADELPKLQVGFIDFVCTFVYKEFSRFHEEITPMLDGITNNRKEWKALADEYETKMKGLEEEKQKQQAANQAAAGSQHGGKQPGGGPASKSCCVQ'
    # Remove whitespaces
    query_sequence = "".join(query_sequence.split())

    # Job name
    str_mode = args.str_mode
    jobname = args.jobname
    template_mode = args.template_mode
    # Remove whitespaces and special characters
    basejobname = "".join(jobname.split())
    basejobname = re.sub(r'\W+', '', basejobname)
    jobname = add_hash(basejobname, query_sequence)

    set_working_directory(args.output_path)

    while os.path.isfile(f"{jobname}.csv"):
        jobname = add_hash(basejobname, ''.join(random.sample(query_sequence, len(query_sequence))))

    with open(f"{jobname}.csv", "w") as text_file:
        text_file.write(f"id,sequence\n{jobname},{query_sequence}")

    queries_path = f"{jobname}.csv"

    # Number of models to use
    num_relax = int(args.num_relax)
    use_amber = num_relax > 0

    # Template mode
    # template_mode = args.template_mode
    #This option is only active for Vanilla AF mode.
    #"none" = no template information is used, "pdb70" = detect templates in pdb70, "custom" - upload and search own templates (PDB or mmCIF format, see [notes below](#custom_templates))
    #TODO: add functions for Vanilla AF
    # custom_template_path = f"{jobname}_template"

    if str_mode == "strategy 1" or str_mode == "strategy 2":
        
        # if not os.path.exists(custom_template_path):
        #     os.mkdir(custom_template_path)
            
        use_templates = True
            
        #     # Assume daq_file is provided as a local path by the user
            
        #     try:
        #         os.rename(daq_file, f"{jobname}_template/1tmp.cif")
        #     except FileNotFoundError:
        #         print("Could not find daq_file to rename.")
        #         return False

        #     return True
            
        template_mode = "custom"

    elif template_mode == "pdb70":
        use_templates = True
        # custom_template_path = None

    elif template_mode == "custom":
        # if not os.path.exists(custom_template_path):
        #     os.mkdir(custom_template_path)
        
        # uploaded_file_path = input("Please enter the local path of the file you wish to upload: ")
        
        use_templates = True
        # filename = os.path.basename(args.)
        # os.rename(uploaded_file_path, f"{jobname}_template/{filename}")

    else:
        custom_template_path = None
        use_templates = False


##MSA part##
def trim_a3m(a3m,daq,good):

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

def ReadA3M(filename):
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

def ReadDAQ(filename,cutoff,dist_cut):
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
    print(daq)
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
    print('HighDAQ',goodpos)

    return daq,goodpos

def save_a3m(file,a3m):
    lines=''
    for name,seq in a3m:
        lines = lines + name+'\n'
        lines = lines + seq+'\n'
    with open(file,'w') as out:
        out.write(lines)

#Use MSA


def msa(args):
    global daq_file, daq_msa, template_mode
    str_mode = args.str_mode
    # pdb_input_path = args.pdb_input_path
    if args.str_mode == "strategy 2":
        # print('Please upload MSA file (a3m format)')

        # Assume cust_msa_file and pdb_input_path are provided as local paths by the user
        cust_msa_file = args.cust_msa_file
        # pdb_input_path = args.pdb_input_path

        if not os.path.isfile(args.pdb_input_path):
            print('Cannot find DAQ-score output file!!')
            exit(1)  # Or handle this case differently
        
        # print(f'User uploaded MSA file at {cust_msa_file}')
        
        try:
            a3m = ReadA3M(args.cust_msa_path)
            daq, good = ReadDAQ(args.pdb_input_path, 0.0, 0.0)
        except FileNotFoundError:
            print("MSA file or DAQ-score output file not found.")
            return False
        
        new_a3m = trim_a3m(a3m, daq, good)
        
        filename = os.path.join(args.output_path, 'trimmed_msa.a3m')
        save_a3m(filename, new_a3m)
        daq_msa = filename
        return daq_msa
    return



# STEP-2: 
def msa_settings(args):
    global daq_msa,query_sequence,queries_path,jobname
    str_mode = args.str_mode
    msa_mode = "mmseqs2_uniref_env"
    pair_mode = "unpaired_paired"
    jobname = args.jobname
    set_working_directory(args.output_path)

    if str_mode == "strategy 2":
        msa_mode = "custom"
        a3m_file = f"{jobname}.custom.a3m"
        if not os.path.isfile(a3m_file):
            custom_msa = daq_msa
            header = 0
            # import fileinput
            for line in fileinput.FileInput(custom_msa,inplace=1):
                if line.startswith(">"):
                    header = header + 1
                if not line.rstrip():
                    continue
                if line.startswith(">") == False and header == 1:
                    query_sequence = line.rstrip()
                print(line, end='')

            os.rename(custom_msa, a3m_file)
            queries_path=a3m_file
            print(f"moving {custom_msa} to {a3m_file}")
    elif msa_mode.startswith("mmseqs2"):
        a3m_file = f"{jobname}.a3m"
    elif msa_mode == "custom":
        a3m_file = f"{jobname}.custom.a3m"
        if not os.path.isfile(a3m_file):
            custom_msa = input("Enter path to the custom MSA file: ")
            header = 0
            for line in fileinput.FileInput(custom_msa, inplace=1):
                if line.startswith(">"):
                    header += 1
                if not line.rstrip():
                    continue
                if not line.startswith(">") and header == 1:
                    query_sequence = line.rstrip()
                print(line, end='')
            os.rename(custom_msa, a3m_file)
            queries_path = a3m_file
            print(f"moving {custom_msa} to {a3m_file}")
    else:
        a3m_file = f"{jobname}.single_sequence.a3m"
        with open(a3m_file, "w") as text_file:
            # query_sequence = input("Enter query sequence: ")
            text_file.write(">1\n%s" % query_sequence)





# Advanced settings
def advanced_setting(args):
    global model_type,num_recycles,recycle_early_stop_tolerance,pairing_strategy,max_msa,num_seeds,use_dropout,save_all,save_recycles,save_to_google_drive,dpi
    model_type = "auto"
    num_recycles = "1"
    recycle_early_stop_tolerance = "auto"
    pairing_strategy = "greedy"

    max_msa = "auto" #@param ["auto", "512:1024", "256:512", "64:128", "32:64", "16:32"]
    num_seeds = 1 #@param [1,2,4,8,16] {type:"raw"}
    use_dropout = False #@param {type:"boolean"}

    num_recycles = None if num_recycles == "auto" else int(num_recycles)
    recycle_early_stop_tolerance = None if recycle_early_stop_tolerance == "auto" else float(recycle_early_stop_tolerance)
    if max_msa == "auto": max_msa = None

    #@markdown #### Save settings
    save_all = False #@param {type:"boolean"}
    save_recycles = False #@param {type:"boolean"}
    save_to_google_drive = False #@param {type:"boolean"}
    #@markdown -  if the save_to_google_drive option was selected, the result zip will be uploaded to your Google Drive
    dpi = 200 #@param {type:"integer"}
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

#@title Run Prediction

def input_features_callback(input_features):
    global display_images
    if display_images:
        plot_msa_v2(input_features)
        plt.show()
        plt.close()

def prediction_callback(protein_obj, length,prediction_result, input_features, mode):
    global display_images
    model_name, relaxed = mode

    if not relaxed:
        if display_images:
            fig = plot_protein(protein_obj, Ls=length, dpi=150)
            plt.show()
            plt.close()

def prediction(args):
    global daq_file, daq_msa, template_mode, query_sequence,use_templates,use_amber,python_version,display_images,jobname,queries_path,pair_mode,msa_mode,num_relax,custom_template_path
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)
    display_images = True #@param {type:"boolean"}
    

    try:
        K80_chk = os.popen('nvidia-smi | grep "Tesla K80" | wc -l').read()
    except:
        K80_chk = "0"
    pass
    if "1" in K80_chk:
        print("WARNING: found GPU Tesla K80: limited to total length < 1000")
    if "TF_FORCE_UNIFIED_MEMORY" in os.environ:
        del os.environ["TF_FORCE_UNIFIED_MEMORY"]
    if "XLA_PYTHON_CLIENT_MEM_FRACTION" in os.environ:
        del os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"]


    # For some reason we need that to get pdbfixer to import
    if use_amber and f"/usr/local/lib/python{python_version}/site-packages/" not in sys.path:
        sys.path.insert(0, f"/usr/local/lib/python{python_version}/site-packages/")



    result_dir = jobname
    log_filename = os.path.join(jobname,"log.txt")
    if not os.path.isfile(log_filename) or 'logging_setup' not in globals():
        setup_logging(Path(log_filename))
    logging_setup = True

    queries, is_complex = get_queries(queries_path)
    model_type = set_model_type(is_complex, model_type)

    if "multimer" in model_type and max_msa is not None:
        use_cluster_profile = False
    else:
        use_cluster_profile = True

    download_alphafold_params(model_type, Path("."))
    results = run(
        queries=queries,
        result_dir=result_dir,
        use_templates=use_templates,
        custom_template_path=custom_template_path,
        num_relax=num_relax,
        msa_mode=msa_mode,
        model_type=model_type,
        num_models=5,
        num_recycles=num_recycles,
        recycle_early_stop_tolerance=recycle_early_stop_tolerance,
        num_seeds=num_seeds,
        use_dropout=use_dropout,
        model_order=[1,2,3,4,5],
        is_complex=is_complex,
        data_dir=Path("."),
        keep_existing_results=False,
        rank_by="auto",
        pair_mode=pair_mode,
        pairing_strategy=pairing_strategy,
        stop_at_score=float(100),
        prediction_callback=prediction_callback,
        dpi=dpi,
        zip_results=False,
        save_all=save_all,
        max_msa=max_msa,
        use_cluster_profile=use_cluster_profile,
        input_features_callback=input_features_callback,
        save_recycles=save_recycles,
    )
    results_zip = f"{jobname}.result.zip"
    os.system(f"zip -r {results_zip} {jobname}")

def get_arguments():
    parser = argparse.ArgumentParser(description='STEP-1: Input Protein Sequence and DAQ result file')

    # Add arguments
    parser.add_argument('--str_mode', type=str, default='strategy 2',
                        help='Select the DAQ-refine strategy. Choices are Vanilla AF2, strategy 1, and strategy 2.',required=True)
    
    parser.add_argument('--query_sequence', type=str, default='',
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

if __name__ == '__main__':
    args = get_arguments()

    # os.environ['PATH'] = f"{os.environ['PATH']}:{maxit_path}"
    os.environ['RCSBROOT']=RCSBROOT
    os.environ['PATH'] += ":RCSBROOT$/bin"
    input_success = get_input(args)
    if not input_success:
        print("Exiting due to error in input.")
        exit(1)
    print("INFO: STEP-1 Input Protein Sequence and DAQ result file started")
    prepare_trimmed_template(args)
    print("Prepare trimmed template finished.")
    daq_msa = msa(args)
    print("MSA finished.(if applicable)")
    print("INFO: STEP-1 Input Protein Sequence and DAQ result file Done")

    print("INFO: STEP-2 Modeling Part started")
    msa_settings(args)
    print("MSA settings finished.")
    advanced_setting(args)
    print("Advanced settings finished.")
    prediction(args)
    print("Prediction finished.")
    print("INFO: STEP-2 Modeling Part Done")



    print("Selected strategy mode:", args.str_mode)
    print("Input query sequence:", args.query_sequence)
    print("Job ID:", args.jobname)
    print("Number of models to use:", args.num_relax)
    print("Template mode:", args.template_mode)
    print("Output directory:", args.output_path)
