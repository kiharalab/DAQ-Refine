#!/bin/bash
echo "INFO : DAQ-refine started"

cd "/bio/kihara-web/www/em/emweb-jobscheduler/algorithms/DAQ"
#source /etc/profile.d/modules.sh
module load miniconda38

#export CUDA_DEVICE_ORDER="PCI_BUS_ID"
#export CUDA_VISIBLE_DEVICES=1
#Inputs
map=$1
structure=$2
output_dir=$3
eval "$(conda shell.bash hook)"
# Acitivate the conda enviroment
module load cryoread

python3 main.py --mode=0 -F=$map -P=$structure --output=$output_dir --window 9 --stride 2 --batch_size=64
python3 writejobyml.py $output_dir

echo "INFO : DAQ-refine Done"