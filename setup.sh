#!/bin/bash

echo "Installation and setup start."

# Set Conda environment path
ENV_PATH="/bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine"

# Create Conda environment at the desired location
module load miniconda38
conda create --prefix $ENV_PATH python=3.8 -y
conda activate $ENV_PATH

# Make sure your conda environment's bin is in PATH
export PATH=$ENV_PATH/bin:$PATH

# Step 1: Download and extract MAXIT
wget https://sw-tools.rcsb.org/apps/MAXIT/maxit-v11.100-prod-src.tar.gz
tar -xzf maxit-v11.100-prod-src.tar.gz

# Step 2: Set Environment
RCSBROOT=$(pwd)/maxit-v11.100-prod-src
export RCSBROOT
export PATH=$PATH:$RCSBROOT/bin
cd $RCSBROOT

# Fetch from kiharalab (based on your Python code)
curl -s 'https://kiharalab.org/emsuites/daq_refine_count.php?pwd=daq_dklab' > /dev/null

# Step 3: Install Dependencies using Conda
conda install -c anaconda bison flex -y
conda install -c conda-forge csh -y
# bash-completion might not be available via conda. You may need to install it globally or find an alternative.

# Step 4: Compile MAXIT
make binary
csh binary.csh
./bin/DictToSdb -ddlFile ./data/ascii/mmcif_ddl.dic -dictFile ./data/ascii/mmcif_pdbx.dic -dictSdbFile mmcif_pdbx.sdb
mv mmcif_pdbx.sdb ./data/binary
rm -f ./bin/DictToSdb ./bin/cif2bin ./bin/connect_main
if [ -e ./mmcif_pdbx.dic-parser.log ]; then
    rm -rf ./mmcif_pdbx.dic-parser.log
fi

# Step 5: Set maxit path
maxit_path=$(pwd)/bin/maxit
export maxit=$maxit_path
cd ..

# Step 6: Install biopython
git clone https://github.com/biopython/biopython
cd biopython
pip install .
cd ..

echo "Installation and setup complete."
