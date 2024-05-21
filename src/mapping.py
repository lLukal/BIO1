from  minimizers import find_minimizers, reverse_complement
from collections import defaultdict

def find_top_frequent_minimizers(minimizer_counts,f):
    sorted_minimizers = sorted(minimizer_counts.items(), key=lambda x: x[1], reverse=True)
    threshold_index = int(f * len(sorted_minimizers))

    minimizers_to_remove = {minimizer for minimizer, _ in sorted_minimizers[:threshold_index]}

    return minimizers_to_remove


def create_minimizer_index(sequence,k,w,f):
    minimizer_counts = defaultdict(int)
    minimizer_positions = defaultdict(list)

    minimizers = find_minimizers(sequence,k,w)

    for minimizer, position, strand, k, w in minimizers:
        minimizer_counts[minimizer] += 1
        minimizer_positions[minimizer].append((position,strand))

    minimizers_to_remove = find_top_frequent_minimizers(minimizer_counts,f)

    for minimizer in list(minimizers_to_remove):
        if minimizer in minimizers_to_remove:
            del minimizer_positions[minimizer]

    return minimizer_positions

def find_matches(fragment, k, w, minimizer_index):
  fragment_minimizers = find_minimizers(fragment, k, w)
 
  matches = []
  rev_matches = []
  seen_matches = set()
  seen_rev_matches = set()
  
  for minimizer, pos, strand, k, w in fragment_minimizers:
    if minimizer in minimizer_index:
      for ref_pos, ref_strand in minimizer_index[minimizer]:
        match = (pos, ref_pos, minimizer, strand, ref_strand)
        if match not in seen_matches:
          matches.append(match)
          seen_matches.add(match)
    elif reverse_complement(minimizer) in minimizer_index:
      for ref_pos, ref_strand in minimizer_index[reverse_complement(minimizer)]:
        match = (pos, ref_pos, minimizer, strand, ref_strand)
        if match not in seen_rev_matches:
          rev_matches.append(match)
          seen_rev_matches.add(match)
 
  return matches, rev_matches