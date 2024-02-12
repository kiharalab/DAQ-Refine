#!/bin/bash
set -e  # Exit on any error
# set -x  # Echo all commands
set -o pipefail  # Exit if any command in a pipeline fails

echo "INFO : DAQ-refine started"

echo "INFO : STEP-0 DAQ started"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ" || { echo "Failed to change directory"; exit 1; }
#source /etc/profile.d/modules.sh
module load miniconda38

#export CUDA_DEVICE_ORDER="PCI_BUS_ID"
#export CUDA_VISIBLE_DEVICES=1
#Inputs
resolution=$1
jobname=$2
chain_id=$3
input_dir=$4
output_dir=$5
map=$6
structure=$7
query_sequence=$8

chain_folder="chain_${chain_id}"
echo "INFO: start DAQ-refine for chain ${chain_id}"
eval "$(conda shell.bash hook)" || { echo "Failed to initialize Conda"; exit 1; }
# Acitivate the conda enviroment
module load cryoread || { echo "Failed to load cryoread"; exit 1; }

# Use the specific python3 interpreter from cryoread environment
CRYOREAD_PYTHON="/apps/miniconda38/envs/cryoread/bin/python3"

# echo $@
$CRYOREAD_PYTHON main.py --mode=0 -F=$map -P=$structure --output="${output_dir}/${chain_folder}" --window 9 --stride 2 --batch_size=64 --server=1  || { echo "main.py failed"; exit 1; }
# $CRYOREAD_PYTHON writejobyml.py $output_dir  || { echo "writejobyml.py failed"; exit 1; }

echo "INFO : STEP-0 DAQ Done"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine" || { echo "Failed to change directory"; exit 1; }
echo "INFO: leave DAQ dir, enter DAQ_refine"
conda deactivate  || { echo "Failed to deactivate Conda environment"; exit 1; }

module load cuda/11.8
conda activate /bio/kihara-web/www/em/emweb-jobscheduler/conda_envs/daq_refine  || { echo "Failed to activate daq_refine environment"; exit 1; }

# which python3
pdb_input_path="${output_dir}/${chain_folder}/daq_score_w9.pdb"
python3 main.py --resolution="$resolution" --jobname="$jobname" --pdb_input_path="$pdb_input_path" --input_path="${input_dir}/${chain_folder}" --output_path="${output_dir}/${chain_folder}" --query_sequence="$query_sequence" || { echo "main.py failed"; exit 1; }

# rerun DAQ
echo "INFO: STEP-4 Computer refined DAQ Started"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ" || { echo "Failed to change directory"; exit 1; }

daqrefined_daq_dir="${output_dir}/${chain_folder}/DAQ"
# daqrefined_structure="${output_dir}/DAQ/input.pdb"
conda deactivate || { echo "Failed to deactivate Conda environment"; exit 1; }
module load cryoread || { echo "Failed to load cryoread"; exit 1; }

# $CRYOREAD_PYTHON main.py --mode=0 -F=$map -P=$daqrefined_structure --output=$daqrefined_output_dir --window 9 --stride 2 --batch_size=512  || { echo "main.py failed"; exit 1; }

#  get the directories end with "s1", "s2", "af2" 
directories=$(find $daqrefined_daq_dir -type d \( -name "*s1" -o -name "*s2" -o -name "*af2" \))

for dir in $directories; do
    daqrefined_structure="$dir/input_relax_0001.pdb" 
    daqrefined_output_dir="$dir"

    $CRYOREAD_PYTHON main.py --mode=0 -F=$map -P=$daqrefined_structure --output=$daqrefined_output_dir --window 9 --stride 2 --batch_size=512 || { echo "main.py failed in $dir"; exit 1; }
done

echo "INFO: STEP-4 Computer refined DAQ Done"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine" || { echo "Failed to change directory"; exit 1; }
echo "INFO: leave DAQ dir, enter DAQ_refine"
yml_dir="${output_dir}/${chain_folder}"
$CRYOREAD_PYTHON writejobyml.py $yml_dir $chain_id  || { echo "writejobyml.py failed"; exit 1; }

