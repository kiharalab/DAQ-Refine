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
strategy="$1 $2"
jobname=$3
pdb_input_path="${5}/daq_score_w9.pdb"
input_dir=$4
output_dir=$5
map=$6
structure=$7
query_sequence=$8
# msa_file=$9

# echo $query_sequence

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

# which python3

python3 main.py --str_mode="$strategy" --jobname="$jobname" --pdb_input_path="$pdb_input_path" --input_path="$input_dir" --output_path="$output_dir" --query_sequence="$query_sequence" || { echo "main.py failed"; exit 1; }

# rerun DAQ
echo "INFO: STEP-4 Computer refined DAQ Started"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ" || { echo "Failed to change directory"; exit 1; }

daqrefined_daq_dir="${output_dir}/DAQ"
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

echo "INFO: STEP-5 Compare and get the best score DAQ result started"
# 

echo "INFO: STEP-5 Compare and get the best score DAQ result Done"


echo "INFO: STEP-6 Visualize structure quality Started"
cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ-Refine" || { echo "Failed to change directory"; exit 1; }
echo "INFO: leave DAQ dir, enter DAQ_refine"

python3 reverse.py $output_dir  || { echo "reverse.py failed"; exit 1; }
echo "INFO: STEP-6 Visualize structure quality Done"


echo "INFO: STEP-7 write job yml Started"

# which python3
$CRYOREAD_PYTHON writejobyml.py $daqrefined_output_dir  || { echo "writejobyml.py failed"; exit 1; }
echo "INFO: STEP-7 write job yml Done"

echo $output_dir

echo "==================================================DAQ-Refine Finished=================================================="
echo "Results stored in: "
echo "$daqrefined_output_dir"
echo "==================================================DAQ-Refine Finished=================================================="

