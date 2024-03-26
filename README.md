# DAQ-Refine


<img src="https://user-images.githubusercontent.com/50850224/184964587-79a4e08d-4edd-4ef8-b69b-dfa8fe3b4804.png" align="left" style="height:240px">

# DAQ-refine: Protein Structure refinement by DAQ-score and ColabFold

## Reference:    
[1][Terashi, G., Wang, X. et al.. Residue-wise local quality estimation for protein models from cryo-EM maps. Nat Methods 19.9 (2022). https://doi.org/10.1038/s41592-022-01574-4](https://www.nature.com/articles/s41592-022-01574-4)   
[2][Terashi, G., Wang, X., Kihara D.. Protein model refinement for cryo-EM maps using AlphaFold2 and the DAQ
score. Acta Cryst. D79 (2022). https://doi.org/10.1107/S2059798322011676](https://doi.org/10.1107/S2059798322011676)   

## Colab instructions
This notebook contains a modified [ColabFold notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) and our tools.
To identify low quality regions in the protein structyure,
Please use [DAQ-score colab](https://colab.research.google.com/drive/1Q-Dj42QjVO8TCOLXMQBJlvm1zInxPkOu?usp=sharing)

Please click "Open In Colab" DAQ-Refine.ipynb and following the instructions to run DAQ-Refine in Colab.

## Local version instructions

### Local Environment Setup with Conda

This guide will walk you through creating a Conda environment for your project and installing the necessary Python and other dependencies. Note that this project is based on DAQ score, so make sure that your have [DAQ-score](https://github.com/kiharalab/DAQ) installed locally.

### Step 1: Create a Conda Environment

First, create a new Conda environment with Python:

```bash
conda create -n daq_refine python=3.9
```

Activate the environment:
```bash
conda activate daq_refine
```

### Step 2: Install Python Dependencies
Install the required Python packages in your Conda environment:

```bash
pip install -r requirements.txt
```
Note: Adjust the package versions according to your project's requirements.

### Step 3: Download and Set Up Maxit
First, enter the DAQ-Refine repository you have installed and change to local branch.
```bash
cd /your/path/to/DAQ-Refine
git checkout local
```

Second, download the Maxit and compile it.
```bash
wget https://sw-tools.rcsb.org/apps/MAXIT/maxit-v11.100-prod-src.tar.gz
tar -xzf maxit-v11.100-prod-src.tar.gz
cd maxit-v11.100-prod-src
make
./bin/DictToSdb -ddlFile ./data/ascii/mmcif_ddl.dic -dictFile ./data/ascii/mmcif_pdbx.dic -dictSdbFile mmcif_pdbx.sdb
mv mmcif_pdbx.sdb ./data/binary
rm -f ./bin/DictToSdb ./bin/cif2bin ./bin/connect_main
if [ -e ./mmcif_pdbx.dic-parser.log ]; then
    rm -rf ./mmcif_pdbx.dic-parser.log;
fi

```
Note: You may need to install additional dependencies like gcc, make, bison, flex, and csh using your system's package manager (e.g., apt for Ubuntu/Debian or yum for CentOS/RedHat).

### Step 4: install Biopython
```bash
cd /your/path/to/install/biopython
git clone https://github.com/biopython/biopython
cd biopython
pip install .
```

### Step 5: install Alphafold/MSA dependencies
```bash
pip install -q --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/kiharalab/ColabFold'
pip install --upgrade dm-haiku
ln -s /you/path/to/python/dist-packages/colabfold colabfold
ln -s /you/path/to/python/dist-packages/alphafold alphafold
sed -i 's/weights = jax.nn.softmax(logits)/logits=jnp.clip(logits,-1e8,1e8);weights=jax.nn.softmax(logits)/g' alphafold/model/modules.py
conda config --set auto_update_conda false
conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 openmm=7.7.0 python='3.9' pdbfixer
conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 python='3.9'
conda install -y -c conda-forge openmm=7.7.0 python='3.9' pdbfixer
```

### Step 6: install MMalign
```bash
cd /your/path/to/DAQ-Refine
wget https://zhanggroup.org/MM-align/bin/module/MMalign.cpp
g++ -static -O3 -ffast-math -o MMalign MMalign.cpp
```

### Step 7: install Rosetta Relaxation
We will use Rosetta Relaxation in the DAQ-Refine final part, so please refer to the [Rosetta](https://www.rosettacommons.org/software/license-and-download) for furture installation


## Usage
### 1. Command parameters
```bash
usage: python3 main.py [-h] [--log_folder_path=LOG_FOLDER_PATH] [--ip_folder_path=IP_FOLDER_PATH] [--op_folder_path=OP_FOLDER_PATH] [--root_run_dir=ROOT_RUN_DIR] [--resolution=RESOLUTION] [--job_id=JOB_ID] [--input_map=INPUT_MAP_PATH] [--pdb_file_path=PDB_FILE_PATH] [--pdb_name=PDB_NAME] [--fasta_file_path=FASTA_FILE_PATH] [--align_strategy=ALIGN_STRATEGY("Manual alignment" or "Smith Waterman")]

required arguments:
  -h, --help               show this help message and exit
  --log_folder_path        log files folder path
  --ip_folder_path         folder path for sanitizing all input files
  --op_folder_path         folder path for output files
  --root_run_dir           root dir for both DAQ git repo and DAQ-Refine git repo
  --input_map              input map path, use the .mrc format, default: input.mrc
  --pdb_file_path          input PDB file path, default: input.pdb
  --fasta_file_path        input fasta file path, default: input.fasta
```






