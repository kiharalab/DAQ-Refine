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
map_path=output_folder+"/input_resize.mrc"
with mrcfile.open(map_path,permissive=True) as mrc:
    data=mrc.data
data=data[data>0]
sort_data=np.sort(data)
contour=float(sort_data[int(0.05*len(sort_data))])
with open(f'{output_folder}/job.yml', 'w') as outfile:
    outfile.write("pdbfiles: daq_score_w9.pdb\n")
    outfile.write("contour: %f\n"%contour)
    outfile.write("strings: DAQ-refine is a protocol using <a href='https://em.kiharalab.org/algorithm/daqscore'>DAQ score</a> to evaluate protein models from cryo-EM maps and employs a modified AlphaFold 2 for refining regions with potential errors. <br> If you encounter any questions for the scored     structure, feel free to email dkihara@purdue.edu, gterashi@purdue.edu and wang3702@purdue.edu.\n")
    #yaml.dump(data, outfile, default_flow_style=False)
