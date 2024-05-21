import argparse
import os
import sys
from Bio import SeqIO

debug = True
def log(message, level=0):
  global debug
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
    log(f"Analyzing file {file_path}...")
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