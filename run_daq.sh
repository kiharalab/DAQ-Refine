#!/bin/bash
#Inputs
map=$1
structure=$2
output_dir=$3
emweb_path=$4
echo "DAQ Started"

cd "$emweb_path/DAQ"
#source /etc/profile.d/modules.sh
# module load miniconda38
#
#export CUDA_DEVICE_ORDER="PCI_BUS_ID"
#export CUDA_VISIBLE_DEVICES=1
eval "$(conda shell.bash hook)"
# Acitivate the conda enviroment
# module load cryoread
which python3
python3 main.py --mode=0 -F=$map -P=$structure --output=$output_dir --window 9 --stride 2 --batch_size=64 --server 1
python3 $emweb_path/DAQ-Refine/utils/gen_score.py $output_dir