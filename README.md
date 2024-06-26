# DAQ-Refine

<img src="https://user-images.githubusercontent.com/50850224/184964587-79a4e08d-4edd-4ef8-b69b-dfa8fe3b4804.png" align="left" width="240">

## DAQ-refine: Refinement of Protein Structures Utilizing DAQ-score and ColabFold

DAQ-Refine utilizes DAQ-score and ColabFold to refine protein structures from cryo-EM maps, aiming to improve prediction accuracy and reliability through detailed residue-wise quality assessment.

## Table of Contents
- [Getting Started with DAQ-Refine](#getting-started-with-daq-refine)
- [DAQ-Refine Comprehensive Guide](#daq-refine-comprehensive-guide)
  - [Source Code Instructions](#source-code-instructions)
  - [Colab Instructions](#colab-instructions)
- [References](#references)


## Getting Started with DAQ-Refine

DAQ-Refine is an advanced protocol designed to evaluate and refine protein models derived from cryo-electron microscopy (cryo-EM) maps. Utilizing a modified AlphaFold2 approach, DAQ-Refine identifies and corrects potential inaccuracies in specified regions of protein models.

### Available Usage Options

Depending on your specific needs and familiarity with DAQ-Refine, we offer three distinct approaches to access and utilize this powerful tool:

### 1. **EM Server:**
Ideal for those new to DAQ-Refine or seeking a straightforward method to obtain results quickly. Visit our [EM Server website](https://em.kiharalab.org/algorithm/daq-refine) for easy access and guidance on submitting your protein models for evaluation and refinement.

### 2. **Source Code:**
For users with unique requirements or those interested in developing a customized version of DAQ-Refine, we provide direct access to our source code. You can download and install it from our [GitHub repository](https://github.com/kiharalab/DAQ-Refine/tree/local). This option allows for extensive customization and integration into existing workflows.

### 3. **Google Colab:**
We also offer a Google Colab notebook for users who prefer an interactive, web-based platform. Access the DAQ-Refine notebook [here](https://colab.research.google.com/github/kiharalab/DAQ-Refine/blob/main/DAQ_Refine.ipynb), upload your protein model, and execute the notebook cells to start the refinement process and receive your results promptly.

## DAQ-Refine Comprehensive Guide

### **Source Code Instructions**

#### Local Environment Setup with Conda

This section guides you through creating a Conda environment for DAQ-Refine, including the installation of Python and other dependencies. Ensure [DAQ-score](https://github.com/kiharalab/DAQ) is installed locally under the same path with DAQ-Refine for optimal performance.

#### Step 1: Create a Conda Environment

Create and activate a new Conda environment:

```bash
conda create -n daq_refine python=3.9
conda activate daq_refine
```

#### Step 2: Install Python Dependencies
Switch to the DAQ-Refine local branch and install the required Python packages in your Conda environment:

```bash
cd /your/path/to/DAQ-Refine
git checkout local
chmode +x install_dependency.sh
./install_dependency.sh
```
*Note:*

- *Adjust the package versions according to your project's requirements.*
- *This script assumes you're in the DAQ-Refine directory to start with.*
- *Adjust paths and commands as necessary, especially if paths differ on your system or if additional configuration is needed.*
- *You may need to install additional dependencies like pytorch, CUDA, gcc, make, bison, flex, Julia and csh using your system's package manager (e.g., apt for Ubuntu/Debian or yum for CentOS/RedHat).*
- *TensorFlow and CUDA configurations may present compatibility issues within your local environment. For comprehensive guidance and resolution strategies, we recommend consulting the official documentation available on their respective websites.*

#### Step 3: install Rosetta Relaxation
We will use Rosetta Relaxation in the DAQ-Refine final part, so please refer to the [Rosetta](https://www.rosettacommons.org/software/license-and-download) for furture installation.


#### Usage
##### 1. Command parameters
```bash
usage: python3 main.py [-h] [--log_folder_path=LOG_FOLDER_PATH] [--ip_folder_path=IP_FOLDER_PATH] [--op_folder_path=OP_FOLDER_PATH] [--root_run_dir=ROOT_RUN_DIR] [--resolution=RESOLUTION] [--job_id=JOB_ID] [--input_map=INPUT_MAP_PATH] [--pdb_file_path=PDB_FILE_PATH] [--pdb_name=PDB_NAME] [--fasta_file_path=FASTA_FILE_PATH] [--align_strategy=ALIGN_STRATEGY("Manual alignment" or "Smith Waterman")] [--rosetta_pth=ROSETTA_PATH]

required arguments:
  -h, --help               show this help message and exit
  --log_folder_path        log files folder path
  --ip_folder_path         folder path for sanitizing all input files
  --op_folder_path         folder path for output files
  --root_run_dir           root dir for both DAQ git repo and DAQ-Refine git repo
  --input_map              input map path, use the .mrc format, default: input.mrc
  --pdb_file_path          input PDB file path, default: input.pdb
  --fasta_file_path        input fasta file path, default: input.fasta
  --rosetta_path           path for the Rosetta script(/your/path/to/rosetta_bin_system_version_bundle/main)
```

### **Colab Instructions**

[This notebook](https://colab.research.google.com/github/kiharalab/DAQ-Refine/blob/main/DAQ_Refine.ipynb) integrates a modified [ColabFold notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) with our DAQ-Refine tools. To pinpoint low-quality regions in the protein structure, please utilize our [DAQ-score Colab notebook](https://colab.research.google.com/drive/1Q-Dj42QjVO8TCOLXMQBJlvm1zInxPkOu?usp=sharing).

To commence with DAQ-Refine in Colab:

- Click [here](https://colab.research.google.com/github/kiharalab/DAQ-Refine/blob/main/DAQ_Refine.ipynb) to open the Colab notebook.
- Follow the instructions within the notebook to execute DAQ-Refine.


### References:

- [1] Terashi, G., Wang, X., Maddhuri Venkata Subramaniya, S. R., Tesmer, J. J., & Kihara, D. (2022). *Residue-wise local quality estimation for protein models from cryo-EM maps*. Nature Methods, 19(9), 1116-1125. [https://doi.org/10.1038/s41592-022-01574-4](https://www.nature.com/articles/s41592-022-01574-4)

- [2] Terashi, G., Wang, X., & Kihara, D. (2023). *Protein model refinement for cryo-EM maps using AlphaFold2 and the DAQ score*. Acta Crystallographica Section D: Structural Biology, 79(1), 10-21. [https://doi.org/10.1107/S2059798322011676](https://doi.org/10.1107/S2059798322011676)