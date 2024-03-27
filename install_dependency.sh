#!/bin/bash

# Step 1: Create and activate a new Conda environment
echo "Creating and activating the Conda environment..."
conda create -n daq_refine python=3.9 -y
conda activate daq_refine

# Assuming you've navigated to the DAQ-Refine directory before running the script
# Step 2: Install Python Dependencies
echo "Switching to local branch and installing Python dependencies..."
git checkout local
pip install -r requirements.txt

# Step 3: Download and Set Up Maxit
echo "Downloading and setting up Maxit..."
wget https://sw-tools.rcsb.org/apps/MAXIT/maxit-v11.100-prod-src.tar.gz
tar -xzf maxit-v11.100-prod-src.tar.gz
cd maxit-v11.100-prod-src
make
./bin/DictToSdb -ddlFile ./data/ascii/mmcif_ddl.dic -dictFile ./data/ascii/mmcif_pdbx.dic -dictSdbFile mmcif_pdbx.sdb
mv mmcif_pdbx.sdb ../data/binary
rm -f ./bin/DictToSdb ./bin/cif2bin ./bin/connect_main
if [ -e ./mmcif_pdbx.dic-parser.log ]; then
    rm -rf ./mmcif_pdbx.dic-parser.log;
fi
cd ..

# Step 4: Install Biopython
echo "Installing Biopython..."
git clone https://github.com/biopython/biopython
cd biopython
pip install .
cd ..

# Step 5: Install Alphafold/MSA dependencies
echo "Installing Alphafold/MSA dependencies..."
pip install -q --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/kiharalab/ColabFold'
pip install --upgrade dm-haiku
ln -s $(python -c 'import site; print(site.getsitepackages()[0])')/colabfold colabfold
ln -s $(python -c 'import site; print(site.getsitepackages()[0])')/alphafold alphafold
sed -i 's/weights = jax.nn.softmax(logits)/logits=jnp.clip(logits,-1e8,1e8);weights=jax.nn.softmax(logits)/g' $(python -c 'import site; print(site.getsitepackages()[0])')/alphafold/model/modules.py
conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 openmm=7.7.0 python='3.9' pdbfixer
conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 python='3.9'

# Step 6: Install MMalign
echo "Downloading and compiling MMalign..."
wget https://zhanggroup.org/MM-align/bin/module/MMalign.cpp
g++ -static -O3 -ffast-math -o MMalign MMalign.cpp

echo "Installation completed."
