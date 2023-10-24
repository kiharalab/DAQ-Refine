#import yaml
import numpy as np
import os
import sys
import glob


# Get the highest score directory
print("INFO: STEP-5 Compare and get the best score DAQ result started")
output_folder = os.path.join(sys.argv[1],"DAQ")
job_folder = sys.argv[1]

def find_highest_score_directory(path):
    # Initialize the highest score and the corresponding directory
    highest_score = 0
    highest_score_directory = ""

    # Get the current path
    current_path = path

    # Iterate through each directory with the specified extensions
    for extension in ['s1', 's2', 'af2']:
        for directory in glob.glob(os.path.join(current_path, f'*_{extension}')):
            # Construct the full path of the file
            file_path = os.path.join(directory, 'daq_raw_score.pdb')

            # Check if the file exists
            if os.path.isfile(file_path):
                with open(file_path, 'r') as file:
                    # Read the last line of the file
                    last_line = file.readlines()[-1]
                    # Parse the last line to get the AA score
                    score = float(last_line.split('AA:')[1].strip())
                    # If this score is higher than the current highest score, update the highest score and the corresponding directory
                    if score > highest_score:
                        highest_score = score
                        highest_score_directory = directory

    return highest_score_directory

# Execute the function and print the result
highest_directory = find_highest_score_directory(output_folder)
print(f'The directory with the highest AA score is: {highest_directory}')

print("INFO: STEP-5 Compare and get the best score DAQ result Done")


# reverse the pdb files
print("INFO: STEP-6 Visualize structure quality Started")
def reverse_pdb(filename, new_file_name):
    with open(filename, "r") as rfile:
        with open(new_file_name, 'w') as wfile:
            for l in rfile:
                if l.startswith('ATOM'):
                    sco = float(l[61:67])
                    line = l[:60] + "%6.2f" % (sco) + "\n"
                    wfile.write(line)

def visualize_structure_quality_3d(output_path):
    refined_result_path = output_path
    final_pdb_path = os.path.join(refined_result_path, "daq_score_w9_reverse.pdb")
    reverse_pdb(os.path.join(refined_result_path, "daq_score_w9.pdb"), final_pdb_path)


def get_score(filename):
    p = {}
    with open(filename) as result:
        for l in result:
            if l.startswith('ATOM') and l[12:16].replace(" ", "") == 'CA':
                l = l.strip("\n")
                split_result = l.split()
                resn = int(float(l[22:26]))
                p[resn] = float(split_result[-1])
    return p

def read_chain_set(filename):
    chain_set = set()
    with open(filename) as result:
        for l in result:
            if l.startswith('ATOM'):
                chain_name = l[21]
                chain_set.add(chain_name)
    return chain_set

def visualize_structure_quality_2d(output_path):
    refined_result_path = os.path.join(output_path,"DAQ")
    output_pdb_path2 = os.path.join(refined_result_path, "daq_score_w9.pdb")
    chain_list = read_chain_set(output_pdb_path2)
    chain_list = list(chain_list)
    chain_list.sort()
    for chain_name in chain_list:
        output_pdb_path3 = os.path.join(refined_result_path, "daq_score_w9_" + str(chain_name) + ".pdb")
        final_pdb_path3 = os.path.join(refined_result_path, "daq_score_w9_" + str(chain_name) + "_reverse.pdb")
        reverse_pdb(output_pdb_path3, final_pdb_path3)

try:
    visualize_structure_quality_3d(highest_directory)
    visualize_structure_quality_2d(highest_directory)
    print("Visualize structure quality finished.")
except Exception as e:
    print(f"Error in visualize_structure_quality(): {e}")
    exit(1)
print( "INFO: STEP-6 Visualize structure quality Done")

# writejobyml
print("INFO: STEP-7 write job yml Started")
data = {
    "pdbfiles":
  "daq_score_w9_reverse.pdb"
}
# parent_folder = os.path.dirname(output_folder)
check_file=highest_directory+"/daq_score_w9_reverse.pdb"
if os.path.exists(check_file):
    with open("%s/done.out"%job_folder,'w') as wfile:
        wfile.write("DONE\n")
else:
    with open("%s/fail.out"%job_folder,'w') as wfile:
        wfile.write("FAIL\n")
import mrcfile

map_path = os.path.join(job_folder, "input_resize.mrc")
with mrcfile.open(map_path,permissive=True) as mrc:
    data=mrc.data
data=data[data>0]
sort_data=np.sort(data)
contour=float(sort_data[int(0.05*len(sort_data))])
with open(f'{job_folder}/job.yml', 'w') as outfile:
    outfile.write("pdbfiles: %s/daq_score_w9_reverse.pdb\n"%highest_directory)
    outfile.write("contour: %f\n"%contour)
    outfile.write("strings: DAQ-refine is a protocol using <a href='https://em.kiharalab.org/algorithm/daqscore'>DAQ score</a> to evaluate protein models from cryo-EM maps and employs a modified AlphaFold 2 for refining regions with potential errors.<br> The 3D model is colored by <a href='https://www.nature.com/articles/s41592-022-01574-4'>DAQ(AA) score</a> scaled from red (-1.0) to blue (1.0) with a 19 residues sliding window. Here blue indicates good score, while red indicates bad score from DAQ. <br> If you encounter any questions for the scored structure, feel free to email dkihara@purdue.edu, gterashi@purdue.edu and wang3702@purdue.edu.\n")
    #yaml.dump(data, outfile, default_flow_style=False)

print("INFO: STEP-7 write job yml Done")

print("==================================================DAQ-Refine Finished==================================================")
print("Results stored in: ")
print(highest_directory)
print("==================================================DAQ-Refine Finished==================================================")