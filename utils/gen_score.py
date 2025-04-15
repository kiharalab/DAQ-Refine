import sys
import os
output_dir=os.path.abspath(sys.argv[1])
input_pdb=os.path.join(output_dir,"daq_score_w9.pdb")
output_txt=output_dir+"/daq_score.txt"
with open(input_pdb,'r') as result:
    with open(output_txt,'w') as wfile:
        for l in result:
            if l.startswith("ATOM") and l[12:16].replace(" ","")=="CA":
                l=l.strip("\n")
                chain_id=l[21]
                resn=int(float(l[22:26]))
                wfile.write("%s %d %.4f\n"%(chain_id,resn,float(l.split()[-1])))