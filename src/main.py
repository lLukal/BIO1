import alignment

import argparse
import os
import sys
from Bio import SeqIO
import datetime
import time
import random


debug = False
def log(message, level=0):
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


######################################################################################################################################################
def main():
  current_time = datetime.datetime.now()
  output_filename = f'./output/output_{current_time.strftime("%Y-%m-%d_%H-%M-%S")}.txt'

  global debug
  args = parse_arguments()
  if args.debug:
    debug = True
    log("Debug mode is enabled.")    
  
  # Load reference genome and fragments files
  reference_filepath = args.reference
  fragments_filepath = args.fragments
  reference, reference_filetype = load_file(reference_filepath)
  fragments, fragments_filetype = load_file(fragments_filepath)
  
  # Adjust filetypes
  if reference_filetype in ['fasta', 'fas', 'fa', 'fna', 'ffn', 'faa', 'mpfa', 'frn']:
    reference_filetype = 'fasta'
  if fragments_filetype in ['fasta', 'fas', 'fa', 'fna', 'ffn', 'faa', 'mpfa', 'frn']:
    fragments_filetype = 'fasta'
  if fragments_filetype in ['fastq', 'fq']:
    fragments_filetype = 'fastq'
  
  log("Input files loaded, filetypes (reference, fragments):")
  log(reference_filetype, 1)
  log(fragments_filetype, 1)
  log("Contents of the files:")
  log(f'{reference[0:10]}...', 1)
  log(f'{fragments[0:10]}...', 1)
  
  stats_reference = analyze(reference_filepath, reference_filetype)
  log(stats_reference[:-1], 1)
  write_output(f'STATS:\n(contig_count, min, max, avg, n_50)\n{str(stats_reference[:-1])}', output_filename)
  stats_fragment = analyze(fragments_filepath, fragments_filetype)
  log(stats_fragment[:-1], 1)
  write_output(f'(contig_count, min, max, avg, n_50)\n{str(stats_fragment[:-1])}\n-------------------\n', output_filename)

  
  reference_seq = stats_reference[-1][0].seq
  fragments_seqs = [record.seq for record in stats_fragment[-1]]
  
  # Align two random fragments
  index_1 = random.randint(0, len(fragments_seqs)-1)
  while len(fragments_seqs[index_1]) > 5000:
    index_1 = random.randint(0, len(fragments_seqs)-1)
  index_2 = random.randint(0, len(fragments_seqs)-1)
  while index_2 == index_1 or len(fragments_seqs[index_2]) > 5000:
    index_2 = random.randint(0, len(fragments_seqs)-1)
  log('')
  log('Found fragments to align...')
  
  log('Computing Global Alignment...', 1)
  start_time = time.perf_counter()
  _, global_score = alignment.needleman_wunsch(fragments_seqs[index_1], fragments_seqs[index_2])
  end_time = time.perf_counter()
  global_time = end_time - start_time

  log('Computing Local Alignment...', 1)
  start_time = time.perf_counter()
  _, local_score = alignment.smith_waterman(fragments_seqs[index_1], fragments_seqs[index_2])
  end_time = time.perf_counter()
  local_time = end_time - start_time
  
  log('Computing Semi-Global Alignment...', 1)
  start_time = time.perf_counter()
  _, semi_score = alignment.semi_global(fragments_seqs[index_1], fragments_seqs[index_2])
  end_time = time.perf_counter()
  semi_time = end_time - start_time
  
  write_output(f'''
Alignment score of two random fragments:
Fragment 1:
\t{fragments_seqs[index_1]}
Fragment 2:
\t{fragments_seqs[index_2]}

Global alignment score: {global_score} (time: {global_time:0.4f}s)
Local alignment score: {local_score} (time: {local_time:0.4f}s)
Semi-global alignment score: {semi_score} (time: {semi_time:0.4f}s)
-------------------
''', output_filename)
  log(f'{global_score, local_score, semi_score}', 1)
  log(f'{global_time, local_time, semi_time}', 1)
  

    
if __name__ == '__main__':
  main()