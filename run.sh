#!/bin/bash
set -e  # Exit on any error
set -o pipefail  # Exit if any command in a pipeline fails

echo "INFO : DAQ-refine started"

echo "INFO : STEP-0 DAQ started"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ" || { echo "Failed to change directory"; exit 1; }
#source /etc/profile.d/modules.sh
module load miniconda38

#export CUDA_DEVICE_ORDER="PCI_BUS_ID"
#export CUDA_VISIBLE_DEVICES=1
#Inputs
strategy=$1
jobname=$2
pdb_input_path="${3}/daq_score_w9.pdb"
input_dir=$3
output_dir=$4
map=$5
structure=$6
# output_dir=$3
eval "$(conda shell.bash hook)" || { echo "Failed to initialize Conda"; exit 1; }
# Acitivate the conda enviroment
module load cryoread || { echo "Failed to load cryoread"; exit 1; }

echo "$4"
python3 main.py --mode=0 -F=$map -P=$structure --output=$output_dir --window 9 --stride 2 --batch_size=64  || { echo "main.py failed"; exit 1; }
python3 writejobyml.py $output_dir  || { echo "writejobyml.py failed"; exit 1; }

echo "INFO : STEP-0 DAQ-refine Done"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine" || { echo "Failed to change directory"; exit 1; }
echo "INFO: leave DAQ dir, enter DAQ_refine"
conda deactivate  || { echo "Failed to deactivate Conda environment"; exit 1; }
conda activate /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine  || { echo "Failed to activate daq_refine environment"; exit 1; }
echo "INFO: STEP-1 Input Protein Sequence and DAQ result file started"

python3 step_1.py --str_mode=$strategy --jobname=$jobname --pdb_input_path=$pdb_input_path --input_path=$input_dir --output_path=$output_dir  || { echo "step_1.py failed"; exit 1; }


echo "INFO: STEP-1 Input Protein Sequence and DAQ result file Done"


