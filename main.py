import os
import sys
from utils.utils import *
from utils.argparser import argparser

def main(args):
    try:
        # Get input parameters from YAML file
        input_map = args.input_map
        input_map = refactor_path(input_map)
        # input_new_map = os.path.join(args.op_folder_path,"input_resize.mrc")
        
        # prepare sequence file
        fasta_file_path = args.fasta_file_path
        fasta_file_path = copy_file(fasta_file_path,os.path.join(args.op_folder_path,"input.fasta"))
        fasta_file_path = format_seq(fasta_file_path,os.path.join(args.op_folder_path,"input_format.fasta"))

        # os.system("python3 utils/reform.py %s %s"%(input_map,input_new_map))
        # input_map = input_new_map
        pdb_file_path = args.pdb_file_path
        pdb_name = args.pdb_name
        pdb_file_path = copy_file(pdb_file_path,os.path.join(args.op_folder_path,"input.pdb"))
        pdb_file_path = correct_pdb_ids(pdb_file_path)

        resolution = args.resolution


        # Extract sequences from PDB file
        pdb_sequences = extract_sequences_from_pdb(pdb_file_path)

        # Parse the FASTA file to get the sequences
        fasta_sequences,nonempty_sequence_length = parse_fasta_file(fasta_file_path)
        
        # if fasta_sequences is False:
        #     print("Error: Sequence is too long, our GPU memory cannot handle it.")
        #     exit()
    except Exception as e:
        print("Error: %s"%e)
        print("Error: Incorrect source files")
        exit()
        

    print("len(fasta_sequences): ",len(fasta_sequences))
    print("len(pdb_sequences): ",len(pdb_sequences))
    # print(fasta_sequences,pdb_sequences)
    if len(fasta_sequences) < len(pdb_sequences):
        print("Error: Mismatch in the number of chains between FASTA and PDB files, fasta length: %d, pdb length: %d"%(len(fasta_sequences),len(pdb_sequences)))
        exit()

    align_strategy = args.align_strategy
    job_id = args.job_id

    # Compare sequences and find the best match
    matches = {}
    try:
        matches = find_best_match(pdb_sequences, fasta_sequences,align_strategy)
    except Exception as e:
        print("Error: %s"%e)
        exit()

    chain_order_file = os.path.join(args.op_folder_path, 'chain_order.txt')

    with open(chain_order_file, 'w') as file:
        for chain_id in matches.keys():
            file.write(chain_id + '\n')

    # Create directories and save the corresponding PDB structures
    with open(pdb_file_path, 'r') as pdb_file:
        pdb_lines = pdb_file.readlines()

    for chain_id, _ in pdb_sequences:
        output_dir = create_directory_for_chain(chain_id,args.log_folder_path)
        output_dir = create_directory_for_chain(chain_id,args.op_folder_path)  
        write_pdb_for_chain(pdb_lines, chain_id, output_dir)

    run_limit = 3600*24*10
    # Run daq score firstly
    daq_1st_file = os.path.join(args.op_folder_path,"daq_score_w9.pdb")
    if check_job_finished(daq_1st_file) == False:
        command_line="bash %s/DAQ-Refine/run_daq.sh "%args.root_run_dir+str(input_map)+" "+str(pdb_file_path)+" "+str(args.op_folder_path) + " " + str(args.root_run_dir)
        os.system(command_line)
        # wait_flag = wait_job(daq_1st_file,run_limit=run_limit)
        # if wait_flag == 2:
        #     print("Error: Time out or daq score failed")
        #     exit()

    chain_num = len(pdb_sequences)
    index = 1
    # begin the work
    for chain_id, (fasta_id, fasta_seq) in matches.items():
        # time.sleep(10)
        if fasta_seq == '' or fasta_seq is None:
            print(f"Warning: Empty sequence found for chain_id {chain_id}")
            index += 1
            continue
        chain_pdb_path = os.path.join('chain_' + chain_id, 'input.pdb')
        chain_pdb_path = os.path.join(args.op_folder_path,chain_pdb_path)
        if chain_pdb_path is None:
            print(f"Warning: No matching pdb file found for chain_id {chain_id}")
            index += 1
            continue  # Skip this chain if no matching pdb file is found
        # Loop to process each chain
        if not os.path.exists(chain_pdb_path):
            print("Dismatch in chain id, fasta_id: %s"%chain_id)
        sequence = preprocess_sequence(fasta_seq)
        structure = chain_pdb_path
        chain_name = f'chain_' + chain_id
        # Begin the processing of the chain
        print("Processing chain %s, index %d/%d"%(chain_id,index,chain_num))
        check_file = os.path.join(args.op_folder_path,chain_name,"done.out")
        fail_file = os.path.join(args.op_folder_path,chain_name,"fail.out")
        log_file = os.path.join(args.log_folder_path,"stdout.log")
        isfinished = check_job_finished(check_file)
        if isfinished:
            print("Chain %s already finished"%chain_id)
            index += 1
            continue
        command_line="bash %s/DAQ-Refine/run_single_chain.sh "%args.root_run_dir+str(resolution)+" "+str(job_id)+" "+str(chain_id)+" "+str(args.ip_folder_path)+" "+str(args.op_folder_path)+" "+str(input_map)+" "+str(structure) + " " + str(sequence) + " " + str(args.root_run_dir)
        # print(command_line)
        os.system(command_line)

        # wait_flag = wait_job(check_file,fail_file,log_file,run_limit=run_limit)
        # print(wait_flag)

        # Copy the log and script files to the output folder
        copy_log_and_script(args.log_folder_path,os.path.join(args.log_folder_path,chain_name))

        # if wait_flag == 2:
        #     print("Error: Time out or chain %s failed"%chain_id)
        #     exit()
        index += 1

    os.system("python3 merge_daqrefine.py %s %s"%(args.op_folder_path,chain_order_file))

if __name__ == '__main__':
    args = argparser()
    main(args)

