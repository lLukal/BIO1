import argparse
import os
import sys
from Bio import SeqIO

debug = False
def log(message, level=0):
  if debug == True:
    indent = ''
    for _ in range(level):
      indent += '\t'
    print(f'{indent}{message}')

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
  if file_type == 'fasta':
    with open(file_path, 'r') as file:
      log(f"Analyzing file {file_path}...")
      contig_lengths = [len(record.seq) for record in SeqIO.parse(file, "fasta")]
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
      
      return len(contig_lengths), min_length, max_length, avg_length, n_50_length


######################################################################################################################################################
def main():
  global debug
  args = parse_arguments()
  if args.debug:
    debug = True
    log("Debug mode is enabled.")    
  
  reference_filepath = args.reference
  fragments_filepath = args.fragments
  # Load reference genome and fragments files
  reference, reference_filetype = load_file(args.reference)
  fragments, fragments_filetype = load_file(args.fragments)
  
  if reference_filetype == 'fna':
    reference_filetype = 'fasta'
  if fragments_filetype == 'fna':
    fragments_filetype = 'fasta'
  
  log("Input files loaded, filetypes (reference, fragments):")
  log(reference_filetype, 1)
  log(fragments_filetype, 1)
  # log("Contents of the files:")
  # log(reference, 0)
  # log(fragments, 0)
  
  stats_reference = analyze(reference_filepath, reference_filetype)
  log(stats_reference, 1)
  stats_fragment = analyze(fragments_filepath, fragments_filetype)
  log(stats_fragment, 1)
  

    
if __name__ == '__main__':
  main()