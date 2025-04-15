import os
import time
import shutil

def timestamp_logging(log_file, message):
    with open(log_file, 'a') as file:
        file.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def fail_job(output_path):
    with open("%s/fail.out"%output_path,'w') as wfile:
        wfile.write("FAIL\n")
    return

def read_fasta_time(input_fasta):
    count=0
    with open(input_fasta,'r') as rfile:
        for line in rfile:
            if line.startswith(">"):
                continue
            count+=len(line)
    return count


def wait_job(check_file, fail_file=None, log_file=None,run_limit=3600 * 24 * 3,check_interval=5,islog = False,timestamp_log=None):
    wait_flag = 0 
    # run_limit = 3600 * 24 * 3  # Set a time limit for the job (3 days)
    start_time = time.time()
    if islog:
        timestamp_logging(timestamp_log, "Start to wait for the job to finish: %s" % check_file)
    # check_interval = 5  # Check every 5 seconds to reduce load on the file system

    while wait_flag == 0:
        try:
            run_time = time.time() - start_time
            if run_time > run_limit:
                print("Time out detected")
                wait_flag = 2
                if islog:
                    timestamp_logging(timestamp_log, "Time out detected: %s" % check_file)
                break  # Exit the loop if time limit is exceeded

            if fail_file and os.path.exists(fail_file):
                print("Failure file detected")
                wait_flag = 2
                if islog:
                    timestamp_logging(timestamp_log, "Failure file detected: %s" % fail_file)
                break  # Exit the loop if failure file is detected

            if log_file and os.path.exists(log_file):
                with open(log_file, 'r') as file:
                    lines = file.readlines()
                    if lines:
                        last_line = lines[-1]
                        if 'fail' in last_line or 'Fail' in last_line:
                            print("Failure detected in log file")
                            wait_flag = 2
                            if islog:
                                timestamp_logging(timestamp_log, "Failure detected in log file: %s" % log_file)
                            break  # Exit the loop if failure is detected in the log file

            if os.path.exists(check_file):
                print("Check file detected, job done.")
                wait_flag = 1
                if islog:
                    timestamp_logging(timestamp_log, "Job done: %s" % check_file)
                break  # Exit the loop if check file is detected

            time.sleep(check_interval)  # Wait for a while before checking again
            if islog:
                timestamp_logging(timestamp_log, "Waiting for the job to finish: %s" % check_file)

        except Exception as e:
            print(f"An error occurred: {e}")
            wait_flag = 2
            if islog:
                timestamp_logging(timestamp_log, "An error occurred: %s" % e)
            break  # Exit the loop if an exception occurs

    return wait_flag


def check_job_finished(check_file):
    if os.path.exists(check_file) and os.path.getsize(check_file) > 0:
        return True
    else:
        return False

import re

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        print (path+" created")
        os.makedirs(path)
        return True
    else:
        print (path+' existed')
        return False

# Created for DAQ-refine preprocessing

def copy_log_and_script(src_path, dest_path):
    """
    Parameters:
    src_path (str): The source directory path.
    dest_path (str): The destination directory path.
    """
    files_to_copy = ['done.out','fail.out']

    # Ensure destination directory exists
    if not os.path.exists(dest_path):
        os.makedirs(dest_path)

    for file_name in files_to_copy:
        src_file_path = os.path.join(src_path, file_name)
        dest_file_path = os.path.join(dest_path, file_name)

        # Check if the source file exists before attempting to copy
        if os.path.exists(src_file_path):
            shutil.copy(src_file_path, dest_file_path)
        else:
            print(f"File not found: {src_file_path}")

# Example usage
# copy_log_and_script('/path/to/source', '/path/to/destination')


def preprocess_sequence(sequence):
    # Remove all non-alphabetic characters
    sequence = re.sub(r'[^a-zA-Z]', '', sequence)
    # Convert all letters to uppercase
    sequence = sequence.upper()
    # Retain only the 20 letters possibly present in protein sequences
    sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence)
    return sequence

def three_to_one(three_letter_code):
    # Convert three-letter code to one-letter code, assuming standard amino acids
    aa_code = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
        'XLE': 'J', 'XAA': 'X', 'UNK': 'X', 'STOP': '*',
    }
    return aa_code.get(three_letter_code, 'X')

def extract_sequences_from_pdb(pdb_path):
    print("Extracting sequences from PDB file...")
    sequences = []
    current_chain = None
    current_sequence = ''
    with open(pdb_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                chain_id = line[21].strip()
                if current_chain != chain_id:
                    if current_sequence:
                        sequences.append((current_chain, current_sequence))
                    current_chain = chain_id
                    current_sequence = ''
                amino_acid = line[17:20].strip()
                current_sequence += three_to_one(amino_acid)
        if current_sequence:
            sequences.append((current_chain, current_sequence))  # Add the last sequence
    return sequences

def parse_fasta_file(fasta_path):
    """
    Parse a FASTA file and return a list of tuples with individual chain IDs and sequences.
    Supports two formats:
    1. Original format with '|': Chain IDs after '|' are split by ','.
    2. New format without '|': Chain IDs are directly split by ','.
    
    :param fasta_path: Path to the FASTA file
    :return: List of tuples (chain ID, sequence) or False if any sequence is longer than the length_limit
    """
    print("Parsing FASTA file...")
    sequences = []
    current_seq = ''
    chain_ids = []
    length_limit = 6000
    nonempty_sequence_length = 0

    with open(fasta_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Save the previous sequence for each chain ID
                if chain_ids:
                    if len(current_seq) > length_limit:
                        return False, 0
                    for chain_id in chain_ids:
                        sequences.append((chain_id.strip(), current_seq))
                        if len(current_seq) > 0:
                            nonempty_sequence_length += 1
                current_seq = ''

                # Determine format and extract chain IDs
                if '|' in line:
                    # Original format with '|'
                    parts = line.split('|')
                    chain_ids = parts[1].strip().split(',') if len(parts) > 1 else ['Unknown']
                else:
                    # New format without '|'
                    chain_ids = line[1:].strip().split(',')

            else:
                current_seq += line.strip()

        # Save the last sequence for each chain ID
        if chain_ids:
            if len(current_seq) > length_limit:
                return False, 0
            for chain_id in chain_ids:
                sequences.append((chain_id.strip(), current_seq))
                if len(current_seq) > 0:
                    nonempty_sequence_length += 1

    return sequences, nonempty_sequence_length


def smith_waterman_alignment(seq1, seq2):
    from Bio.Align import substitution_matrices
    from Bio import Align
    # Create an aligner object
    aligner = Align.PairwiseAligner()

    # Set the aligner to perform local alignment (Smith-Waterman)
    aligner.mode = 'local'
    
    # Load the BLOSUM62 substitution matrix
    blosum62 = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = blosum62

    # Set gap scores: opening and extension
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # Perform the alignment
    alignments = aligner.align(seq1, seq2)

    # not enough similarity
    if not alignments:
        return None

    # Return the best alignment
    best_alignment = alignments[0]
    return best_alignment

def create_directory_for_chain(chain_id,output_folder):
    dir_name = f'chain_{chain_id}'
    dir_name = os.path.join(output_folder,dir_name)
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def write_pdb_for_chain(pdb_lines, chain_id, output_dir):
    chain_lines = [line for line in pdb_lines if line.startswith("ATOM") and line[21].strip() == chain_id]
    with open(f'{output_dir}/input.pdb', 'w') as chain_file:
        chain_file.writelines(chain_lines)

def find_best_match(pdb_sequences, fasta_sequences, strategy=None):
    try:
        from Bio.Seq import Seq
    except ImportError:
        print("Please install the Biopython package to use this algorithm.")
        return {}  
    
    matches = {}
    matched_fasta_indices = set()  # To track the indices of matched FASTA sequences
   

    if strategy == "Manual alignment":
        # Manually match sequences in the order they appear
        for i, (pdb_id, _) in enumerate(pdb_sequences):
            fasta_id, fasta_seq = fasta_sequences[i] 
            matches[pdb_id] = (fasta_id, fasta_seq)

    
    else:
        # using Smith
        score_matrix = [[-float('inf')] * len(pdb_sequences) for _ in range(len(fasta_sequences))]
        # calculate the score matrix
        for i, (_, fasta_seq) in enumerate(fasta_sequences):
            if len(fasta_seq) == 0:  # if the fasta sequence is empty,
                continue
            fasta_seq_obj = Seq(fasta_seq)
            for j, (_, pdb_seq) in enumerate(pdb_sequences):
                pdb_seq_obj = Seq(pdb_seq)
                # calculate the score of the alignment
                alignment = smith_waterman_alignment(pdb_seq_obj, fasta_seq_obj)
                # print(alignment)
                if alignment:
                    score_matrix[i][j] = alignment.score  # Get the score of the alignment

        
        match_candidates = []
        for i, row in enumerate(score_matrix):
            for j, score in enumerate(row):
                if score != -float('inf'):
                    match_candidates.append((score, i, j))
        match_candidates.sort(reverse=True)  # sort the match candidates by score

        used_pdb = set()
        used_fasta = set()

        for score, i, j in match_candidates:
            if j not in used_pdb and i not in used_fasta:
                pdb_id = pdb_sequences[j][0]
                fasta_id, fasta_seq = fasta_sequences[i]
                matches[pdb_id] = (fasta_id, fasta_seq)
                used_pdb.add(j)
                used_fasta.add(i)

    return matches

def check_and_merge_multiple_pdb_files(file_paths, output_file_path):
    """
    Merges multiple PDB files into one. Adds MODEL and ENDMDL tags if they are not present in the original files.

    Parameters:
    file_paths (list): List of paths to the PDB files to be merged.
    output_file_path (str): Path for the output merged PDB file.
    """
    def add_model_tags_if_needed(file_content, model_number):
        if "MODEL" not in file_content and "ENDMDL" not in file_content:
            return f"MODEL        {model_number}\n{file_content}ENDMDL\n"
        else:
            return file_content

    with open(output_file_path, 'w') as output_file:
        for idx, file_path in enumerate(file_paths, start=1):
            with open(file_path, 'r') as file:
                file_content = file.read()
                file_content = add_model_tags_if_needed(file_content, idx)
                output_file.write(file_content)

# Example usage
# check_and_merge_multiple_pdb_files(['path_to_first_pdb_file.pdb', 'path_to_second_pdb_file.pdb', ...], 'path_to_output_merged_file.pdb')

# This function doesn't return anything, it creates a new file with the merged content and adds separators if needed.

def format_chain_name_list():
    # Generate a list of potential chain IDs (single and double-letter codes, digits)
    tmp_chain_list = [chr(i) for i in range(ord('A'), ord('Z') + 1)] + \
                     [chr(i) for i in range(ord('a'), ord('z') + 1)] + \
                     [str(i) for i in range(10)]
    extend_chain_list = [a + b for a in tmp_chain_list for b in tmp_chain_list]
    return tmp_chain_list + extend_chain_list

def correct_pdb_ids(pdb_file_path):
    try:
        with open(pdb_file_path, 'r') as file:
            lines = file.readlines()
    except IOError:
        print("Error: File not accessible or not found.")
        return

    corrected_lines = []
    chain_ids = set()
    temp_ids = format_chain_name_list()
    temp_index = 0
    last_chain_id = None

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
            if line.startswith("TER"):
                try:
                    chain_id = line[21].strip()
                except IndexError:
                    continue
            else:
                chain_id = line[21]

            if chain_id != last_chain_id:
                new_chain_id = None
                last_chain_id = chain_id  # update last_chain_id

                # assign new chain ID if needed
                if chain_id not in temp_ids or chain_id in chain_ids:
                    while temp_ids[temp_index] in chain_ids:
                        temp_index += 1
                        if temp_index >= len(temp_ids):
                            print("Error: Ran out of unique chain IDs.")
                            return
                    new_chain_id = temp_ids[temp_index]
                    temp_index += 1
                    chain_ids.add(new_chain_id)
                    line = line[:21] + new_chain_id + line[22:]
                else:
                    chain_ids.add(chain_id)               

            else:
                if chain_id not in temp_ids or new_chain_id is not None:
                    line = line[:21] + new_chain_id + line[22:]    

        # print(line)
        corrected_lines.append(line)

    new_file_path = pdb_file_path.replace('.pdb', '_corrected.pdb')
    try:
        with open(new_file_path, 'w') as corrected_file:
            corrected_file.writelines(corrected_lines)
        print("PDB file corrected and saved as", new_file_path)
    except IOError:
        print("Error: File could not be written.")
    
    return new_file_path


# Created for DAQ-refine preprocessing

def refactor_path(path):
    new_path="'%s'"%path
    return new_path


# def copy_file(old_path,new_path):
#     shutil.copy(old_path,new_path)
#     return  new_path

def check_allaa(unet_dir):
    label_list = ["ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU", "ASP", "GLU", "LYS", "ARG", "SER", "THR", "TYR",
                        "HIS", "CYS", "ASN", "TRP", "GLN", "GLY"]
    pre_name = "sigmoidAA"
    for k, base_name in enumerate(label_list):
        cur_map_path = os.path.join(unet_dir, pre_name+"_" + str(base_name) + ".mrc")
        if not os.path.exists(cur_map_path):
            return False
    return True


def format_seq(input_path,output_path):
    AA20='ARNDCQEGHILKMFPSTWYV'

    outlines=''
    seq=''
    with open(input_path,'r') as f:
        lines = [l for l in f]
        #print(lines)
        for l in lines:
            l = l.strip()
            #Header lines
            if l.startswith('>'):
                if len(seq)>0:
                    outlines += f'{seq}\n'
                outlines += f'{l}\n'
                seq=''
            else:
                #check 20AA or not
                for aa in l:
                    aa=aa.upper()
                    if aa in AA20:
                        seq+=aa
                        #print(aa)
                    else:
                        print(f'Unknown AA type {aa}')
        if len(seq)>0:
            outlines += f'{seq}\n'
        print(outlines)
    with open(output_path,'w') as out:
        out.write(outlines)
    return output_path
import os

def copy_file(source_path, destination_path):
  """Copies a file to a new path.

  Args:
    source_path: The path to the source file.
    destination_path: The path to the destination file.
  """

  if not os.path.exists(os.path.dirname(destination_path)):
      os.makedirs(os.path.dirname(destination_path), exist_ok=True)

  shutil.copyfile(source_path, destination_path)
  return destination_path



