import os
import sys

def reverse_pdb(filename, new_file_name):
    with open(filename, "r") as rfile:
        with open(new_file_name, 'w') as wfile:
            for l in rfile:
                if l.startswith('ATOM'):
                    sco = float(l[61:67])
                    line = l[:60] + "%6.2f" % (sco) + "\n"
                    wfile.write(line)

def visualize_structure_quality_3d(output_path):
    refined_result_path = os.path.join(output_path,"DAQ")
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

output_path = sys.argv[1]

try:
    visualize_structure_quality_3d(output_path)
    visualize_structure_quality_2d(output_path)
    print("Visualize structure quality finished.")
except Exception as e:
    print(f"Error in visualize_structure_quality(): {e}")
    exit(1)