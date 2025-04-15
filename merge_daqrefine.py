import os
import sys
import yaml
import numpy as np
from utils.utils import check_and_merge_multiple_pdb_files


def merge_pdb_files_in_folders(base_path, output_file, chain_order_file):
    """
    Merge PDB files specified in job.yml files located in folders starting with 'chain_'
    into a single PDB file with modified chain identifiers and separation lines,
    based on the chain order provided in a file.

    Parameters:
    base_path (str): Path to the directory containing the chain folders.
    output_file (str): Path to the output PDB file.
    chain_order_file (str): Path to the file containing the order of chain IDs.

    Returns:
    None
    """
    # Read the chain order from the file
    with open(chain_order_file, 'r') as order_file:
        chain_ids = [line.strip() for line in order_file]

    with open(output_file, 'w') as outfile:
        for chain_id in chain_ids:
            folder = f'chain_{chain_id}'
            folder_path = os.path.join(base_path, folder)

            if os.path.isdir(folder_path):
                job_file_path = os.path.join(folder_path, 'job.yml')
                if os.path.isfile(job_file_path):
                    with open(job_file_path, 'r') as ymlfile:
                        try:
                            cfg = yaml.safe_load(ymlfile)
                            pdb_file_path = cfg['pdbfiles']
                            pdb_file_path = os.path.join(folder_path, pdb_file_path)

                            if pdb_file_path and os.path.isfile(pdb_file_path):
                                with open(pdb_file_path, 'r') as infile:
                                    for line in infile:
                                        if line.startswith("ATOM"):
                                            line = line[:21] + chain_id + line[22:]
                                        outfile.write(line)

                                    # Write separation line
                                    outfile.write("TER\n")
                        except yaml.YAMLError as exc:
                            print(f"Error in configuration file: {exc}")
        # Write the 'END' line at the end of the merged file
        outfile.write("END\n")

base_path = sys.argv[1]  # Replace with the path to your folders
chain_order_file = sys.argv[2] # Replace with the path to your chain order file
output_file = os.path.join(base_path,'merged_refined_protein.pdb')  # The name of the output file
merge_pdb_files_in_folders(base_path, output_file, chain_order_file)
pdb_list = [base_path+"/daq_score_w9.pdb",output_file]
print(pdb_list)
check_and_merge_multiple_pdb_files(pdb_list,base_path+"/merged_protein.pdb")

# writetojobyml
data = {
    "pdbfiles":
  "merged_protein.pdb"
}
# parent_folder = os.path.dirname(output_folder)
check_file=base_path+"/merged_protein.pdb"
if os.path.exists(check_file):
    with open("%s/done.out"%base_path,'w') as wfile:
        wfile.write("DONE\n")
else:
    with open("%s/fail.out"%base_path,'w') as wfile:
        wfile.write("FAIL\n")
import mrcfile

map_path = os.path.join(base_path, "input_resize.mrc")
with mrcfile.open(map_path,permissive=True) as mrc:
    data=mrc.data
data=data[data>0]
sort_data=np.sort(data)
contour=float(sort_data[int(0.05*len(sort_data))])
with open(f'{base_path}/job.yml', 'w') as outfile:
    outfile.write("pdbfiles: merged_protein.pdb\n")
    outfile.write("contour: %f\n"%contour)
    outfile.write("strings: DAQ-refine is a protocol using <a href='https://em.kiharalab.org/algorithm/daqscore'>DAQ score</a> to evaluate protein models from cryo-EM maps and employs a modified AlphaFold 2 for refining regions with potential errors.<br> The 3D model is colored by <a href='https://www.nature.com/articles/s41592-022-01574-4'>DAQ(AA) score</a> scaled from red (-1.0) to blue (1.0) with a 19 residues sliding window. Here blue indicates good score, while red indicates bad score from DAQ. <br>In MODEL1, it shows the original model, all modeled positions are colored by DAQ(AA) score. <br>In MODEL2, it shows the DAQ refined model, all modeled positions are colored by DAQ(AA) score. <br> If you encounter any questions for the scored structure, feel free to email dkihara@purdue.edu, gterashi@purdue.edu and wang3702@purdue.edu.\n")
    #yaml.dump(data, outfile, default_flow_style=False)
