import argparse
import os
import sys
from Bio import SeqIO

def log(message, level=0):
  debug = parse_arguments().debug
  if debug == True:
    indent = ''
    for _ in range(level):
      indent += '\t'
    print(f'{indent}{message}')
    
def write_output(message, output_filename='./output/output_temp.txt'):
  with open(output_filename, 'a') as file:
    file.write(str(message) + '\n')
    
def parse_arguments():
  parser = argparse.ArgumentParser(description="Process reference genome and fragments files.")
  parser.add_argument('--reference', required=True, help="Path to the reference genome file (FASTA format).")
  parser.add_argument('--fragments', required=True, help="Path to the fragments file (FASTA or FASTQ format).")
  parser.add_argument('--debug', action='store_true', help="Enable debug mode.")
  parser.add_argument('--cigar', action='store_true', help="Include the CIGAR string in the output.")
  parser.add_argument('--threads', type=int, help="Number of threads to use.")
  parser.add_argument('--k', type=int, help="Value for parameter k.")
  parser.add_argument('--w', type=int, help="Value for parameter w.")
  parser.add_argument('--f', type=float, help="Value for parameter f.")
  
  return parser.parse_args()

def load_file(file_path):
  if not os.path.isfile(file_path):
    print(f"Error: file '{file_path}' does not exist.", file=sys.stderr)
    sys.exit(1)
  file_extension = file_path.split('.')[-1]
  with open(file_path, 'r') as file:
    return file.read(), file_extension

def analyze(file_path, file_type):
  with open(file_path, 'r') as file:
    log(f"Analyzing file {file_path}...", 1)
    contigs = list(SeqIO.parse(file, file_type))
    contig_lengths = [len(record.seq) for record in contigs]
    max_length = max(contig_lengths)
    min_length = min(contig_lengths)
    avg_length = sum(contig_lengths) / len(contig_lengths)
    
    sorted_lengths = sorted(contig_lengths, reverse=False)
    total_length = sum(sorted_lengths)
    half_total_length = total_length / 2
    cumulative_length = 0
    for length in sorted_lengths:
      cumulative_length += length
      if cumulative_length >= half_total_length:
        n_50_length = length
        break
    
    return len(contig_lengths), min_length, max_length, avg_length, n_50_length, contigs
  
  
###########################################################################################################################3

def calculate_mapping_length(seq1, seq2):
    """Calculate the mapping length (number of bases that are aligned and not gaps)."""
    return sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')

def calculate_number_of_matches(seq1, seq2):
    """Calculate the number of matching bases."""
    return sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')

def write_paf(seq1, seq2, output_file, query_name="query", target_name="target"):
    query_length = len(seq1.replace('-', ''))
    target_length = len(seq2.replace('-', ''))
    mapping_length = calculate_mapping_length(seq1, seq2)
    number_of_matches = calculate_number_of_matches(seq1, seq2)
    mapping_quality = 60  # Placeholder for mapping quality, this usually requires a more complex calculation

    # Find the start and end positions
    query_start = seq1.find(seq1.lstrip('-'))
    query_end = len(seq1.rstrip('-'))
    target_start = seq2.find(seq2.lstrip('-'))
    target_end = len(seq2.rstrip('-'))

    with open(output_file, 'w') as f:
        f.write(f"#query_name\t#query_length\t#query_start\t#query_end\t+\t")
        f.write(f"#target_name\t#target_length\t#target_start\t#target_end\t")
        f.write(f"#number_of_matches\t#mapping_length\t#mapping_quality\n")
        
        f.write(f"{query_name}\t{query_length}\t{query_start}\t{query_end}\t+\t")
        f.write(f"{target_name}\t{target_length}\t{target_start}\t{target_end}\t")
        f.write(f"{number_of_matches}\t{mapping_length}\t{mapping_quality}\n")
