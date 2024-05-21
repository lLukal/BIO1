import mapping
from misc import log, write_output, parse_arguments, load_file, analyze
from visualization import plot_mapped_genome
import datetime
import random

def align_fragment(fragment, reference_seq, minimizer_index, fragment_minimizers, k, w, f):
  log('Finding matches...', 1)
  matches, rev_matches = mapping.find_matches(fragment, k, w, f, minimizer_index, fragment_minimizers)

  log('Finding longest increasing subsequence...', 1)
  lis = mapping.longest_increasing_subsequence(matches)
  log('Finding longest increasing subsequence (complement)...', 1)
  rev_lis = mapping.longest_increasing_subsequence(rev_matches)
  lis = lis if len(lis) > len(rev_lis) else rev_lis
  
  q_begin, q_end = lis[0][0], lis[-1][0]
  t_begin, t_end = lis[0][1], lis[-1][1]

  log('Aligning region...', 1)
  log(f'fragment_begin: {q_begin}, fragment_end: {q_end}, reference_begin: {t_begin}, reference_end: {t_end}', 2)
  log(f'fragment length: {len(fragment)}', 2)
  aligned_seq1, aligned_seq2, alignment_score = mapping.align_region(fragment[q_begin:q_end], reference_seq[t_begin:t_end])
  
  return fragment, 'reference_genome', q_begin, q_end, t_begin, t_end, aligned_seq1, aligned_seq2, alignment_score

def main():
  current_time = datetime.datetime.now()
  output_filename = f'./output/output_{current_time.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
  args = parse_arguments()
  
  k = args.k if args.k else 5
  w = args.w if args.w else 15
  f = args.f if args.f else 0.001
  debug = args.debug
  threads = args.threads
  cigar = args.cigar
  reference_filepath = args.reference
  fragments_filepath = args.fragments

  log(f'--Using: k={k}, w={w}, f={f} debug={debug} threads={threads} cigar={cigar}', 0)
  
  # Load reference genome and fragments files
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
  log(f'Reference: {reference_filetype}, Fragments: {fragments_filetype}', 1)
  
  stats_reference = analyze(reference_filepath, reference_filetype)
  write_output(f'STATS:\n(contig_count, min, max, avg, n_50)\n{str(stats_reference[:-1])}', output_filename)
  stats_fragment = analyze(fragments_filepath, fragments_filetype)
  write_output(f'(contig_count, min, max, avg, n_50)\n{str(stats_fragment[:-1])}\n-------------------\n', output_filename)

  #############################################################################
  
  reference_seq = stats_reference[-1][0].seq
  all_fragments_seqs = [record.seq for record in stats_fragment[-1]]
  
  fragments_seqs = []
  while len(fragments_seqs) < 4:
    random_int = random.randint(0, len(stats_fragment[-1]))
    if len(all_fragments_seqs[random_int]) <= 5000:
      fragments_seqs.append(all_fragments_seqs[random_int])
  
  #############################################################################
  results = []
  
  log('Creating minimizer index for reference...')
  minimizer_index = mapping.create_minimizer_index(reference_seq, k, w, f)
  log('--------------------------------------')
  
  for index, fragment in enumerate(fragments_seqs):
    log(f'Fragment {index+1} of {len(fragments_seqs)}...')
    log('Finding minimizers for fragment...', 1)
    fragment_minimizers = mapping.create_minimizer_index(fragment, k, w, f)
    
    result = align_fragment(fragment, reference_seq, minimizer_index, fragment_minimizers, k, w, f)
    results.append(result)
    log(f'Done!', 1)
  
  # write_output('\n\nPAF:\n', output_filename)
  # with ThreadPoolExecutor(max_workers=4) as executor:
  #   future_to_fragment = {executor.submit(align_fragment, fragment, reference_seq, minimizer_index): fragment for fragment in fragments_seqs}
  #   for future in as_completed(future_to_fragment):
  #     try:
  #       result = future.result()
  #       results.append(result)
  #     except Exception as exc:
  #       log(f'Fragment generated an exception: {exc}')
  
  # Write results to the output file
  for result in results:
    # fragment, ref_name, q_begin, q_end, t_begin, t_end, aligned_seq1, aligned_seq2, alignment_score = result
    # log('Printing PAF...')
    # mapping.print_paf(fragment, ref_name, q_begin, q_end, t_begin, t_end, aligned_seq1, aligned_seq2, alignment_score)
    # write_output(f'{fragment}\t{len(fragment)}\t{q_begin}\t{q_end}\t+\t{ref_name}\t{len(ref_name)}\t{t_begin}\t{t_end}\t{alignment_score}\t60', output_filename)
    pass
  
  plot_mapped_genome(results, len(reference_seq))


if __name__ == '__main__':
  main()