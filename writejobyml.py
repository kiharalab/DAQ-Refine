#import yaml
import numpy as np
import os
import sys

output_folder = sys.argv[1]

data = {
    "pdbfiles":
  "daq_score_w9.pdb"
}
check_file=output_folder+"/daq_score_w9.pdb"
if os.path.exists(check_file):
    with open("%s/done.out"%output_folder,'w') as wfile:
        wfile.write("DONE\n")
else:
    with open("%s/fail.out"%output_folder,'w') as wfile:
        wfile.write("FAIL\n")
import mrcfile
parent_folder = os.path.dirname(output_folder)
map_path = os.path.join(parent_folder, "input_resize.mrc")
with mrcfile.open(map_path,permissive=True) as mrc:
    data=mrc.data
data=data[data>0]
sort_data=np.sort(data)
contour=float(sort_data[int(0.05*len(sort_data))])
with open(f'{output_folder}/job.yml', 'w') as outfile:
    outfile.write("pdbfiles: daq_score_w9.pdb\n")
    outfile.write("contour: %f\n"%contour)
    outfile.write("strings: DAQ-refine is a protocol using <a href='https://em.kiharalab.org/algorithm/daqscore'>DAQ score</a> to evaluate protein models from cryo-EM maps and employs a modified AlphaFold 2 for refining regions with potential errors.<br> The 3D model is colored by <a href='https://www.nature.com/articles/s41592-022-01574-4'>DAQ(AA) score</a> scaled from red (-1.0) to blue (1.0) with a 19 residues sliding window. Here blue indicates good score, while red indicates bad score from DAQ. <br> If you encounter any questions for the scored structure, feel free to email dkihara@purdue.edu, gterashi@purdue.edu and wang3702@purdue.edu.\n")
    #yaml.dump(data, outfile, default_flow_style=False)
