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
strategy="$1 $2"
jobname=$3
pdb_input_path="${5}/daq_score_w9.pdb"
input_dir=$4
output_dir=$5
map=$6
structure=$7
eval "$(conda shell.bash hook)" || { echo "Failed to initialize Conda"; exit 1; }
# Acitivate the conda enviroment
module load cryoread || { echo "Failed to load cryoread"; exit 1; }

# Use the specific python3 interpreter from cryoread environment
CRYOREAD_PYTHON="/apps/miniconda38/envs/cryoread/bin/python3"

# echo $@
$CRYOREAD_PYTHON main.py --mode=0 -F=$map -P=$structure --output=$output_dir --window 9 --stride 2 --batch_size=64  || { echo "main.py failed"; exit 1; }
# $CRYOREAD_PYTHON writejobyml.py $output_dir  || { echo "writejobyml.py failed"; exit 1; }

echo "INFO : STEP-0 DAQ-refine Done"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine" || { echo "Failed to change directory"; exit 1; }
echo "INFO: leave DAQ dir, enter DAQ_refine"
conda deactivate  || { echo "Failed to deactivate Conda environment"; exit 1; }

module load cuda/11.8
conda activate /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine  || { echo "Failed to activate daq_refine environment"; exit 1; }

which python3

python3 main.py --str_mode="$strategy" --jobname="$jobname" --pdb_input_path="$pdb_input_path" --input_path="$input_dir" --output_path="$output_dir"  || { echo "main.py failed"; exit 1; }

daqrefined_output_dir="${output_dir}/DAQ"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine" || { echo "Failed to change directory"; exit 1; }
echo "INFO: leave DAQ dir, enter DAQ_refine"
python3 writejobyml.py $daqrefined_output_dir  || { echo "writejobyml.py failed"; exit 1; }



