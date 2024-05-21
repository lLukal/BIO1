from minimizers import find_minimizers, reverse_complement
from collections import defaultdict, deque
from alignment import needleman_wunsch, smith_waterman, semi_global
from bisect import bisect_left
from misc import log

def find_top_frequent_minimizers(minimizer_counts, f):
  sorted_minimizers = sorted(minimizer_counts.items(), key=lambda x: x[1], reverse=True)
  threshold_index = int(f * len(sorted_minimizers))
  
  minimizers_to_remove = {minimizer for minimizer, _ in sorted_minimizers[:threshold_index]}
  
  return minimizers_to_remove

def create_minimizer_index(sequence, k, w, f):
  minimizer_counts = defaultdict(int)
  minimizer_positions = defaultdict(list)
  
  minimizers = find_minimizers(sequence, k, w)
  for minimizer, pos, strand, k, w in minimizers:
    minimizer_counts[minimizer] += 1
    minimizer_positions[minimizer].append((pos, strand))
  
  minimizers_to_remove = find_top_frequent_minimizers(minimizer_counts, f)
  
  for minimizer in list(minimizer_positions.keys()):
    if minimizer in minimizers_to_remove:
      del minimizer_positions[minimizer]
  
  return minimizer_positions

def find_matches(fragment, k, w, f, minimizer_index, fragment_minimizers = None):
  if fragment_minimizers is None:
    fragment_minimizers = find_minimizers(fragment, k, w)

  matches = []
  rev_matches = []
  # seen_matches = set()
  # seen_rev_matches = set()
  
  for minimizer, appearances in fragment_minimizers.items():
    for appearance in appearances:
      pos, strand = appearance
      if minimizer in minimizer_index:
        for ref_pos, ref_strand in minimizer_index[minimizer]:
          match = (pos, ref_pos, minimizer, strand, ref_strand)
          # if match not in seen_matches:
          matches.append(match)
          # seen_matches.add(match)
      elif reverse_complement(minimizer) in minimizer_index:
        for ref_pos, ref_strand in minimizer_index[reverse_complement(minimizer)]:
          match = (pos, ref_pos, minimizer, strand, ref_strand)
          # if match not in seen_rev_matches:
          rev_matches.append(match)
          # seen_rev_matches.add(match)

  matches = sorted(matches, key=lambda x: (x[0], x[1]))
  rev_matches = sorted(rev_matches, key=lambda x: (x[0], x[1]))
  
  return matches, rev_matches
  
def longest_increasing_subsequence(matches):
  positions = [match[0] for match in matches]

  n = len(positions)
  lis = []
  predecessors = [-1] * n
  lis_indices = []

  for i, pos in enumerate(positions):
    idx = bisect_left([positions[j] for j in lis_indices], pos)
    
    if idx < len(lis_indices):
      lis_indices[idx] = i
    else:
      lis_indices.append(i)
    
    if idx > 0:
      predecessors[i] = lis_indices[idx - 1]
    else:
      predecessors[i] = -1
      
    if i % 100000 == 0:
      log(f"{i} / {n}", 2)
  
  if not lis_indices:
    return []
  
  k = lis_indices[-1]
  lis = []
  while k != -1:
    lis.append(matches[k])
    k = predecessors[k]
  
  lis.reverse()
  return lis

def align_region(fragment, ref_seq):
  # Use Needleman-Wunsch alignment algorithm on the specified region
  _, aligned_seq1, aligned_seq2, alignment_score = needleman_wunsch(fragment, ref_seq)
  return aligned_seq1, aligned_seq2, alignment_score

def print_paf(fragment, ref_name, q_begin, q_end, t_begin, t_end, aligned_seq1, aligned_seq2, alignment_score):
#   print(f"{fragment}\t{len(fragment)}\t{q_begin}\t{q_end}\t+\t{ref_name}\t{len(ref_name)}\t{t_begin}\t{t_end}\t{alignment_score}\t60")
  with open('./output/output.paf', 'a') as file:
    file.write(f"{fragment}\t{len(fragment)}\t{q_begin}\t{q_end}\t+\t{ref_name}\t{len(ref_name)}\t{t_begin}\t{t_end}\t{alignment_score}\t60\n")
    file.write(f"{aligned_seq1}\n")
    file.write(f"{aligned_seq2}\n")
    file.write(f"\n")
