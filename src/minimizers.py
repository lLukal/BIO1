import multiprocessing as mp

def reverse_complement(seq):
  complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
  return ''.join(complement[base] for base in reversed(seq))

def find_window_minimizer(args):
  sequence, k, w, i, offset = args
  window = sequence[i:i+w+k-1]
  min_k_mer = window[:k]
  min_pos = 0
  min_strand = "original"
  
  for j in range(1, w - k + 1):
    k_mer = window[j:j+k]
    k_mer_rc = reverse_complement(k_mer)
    
    if k_mer < min_k_mer:
      min_k_mer = k_mer
      min_pos = j
      min_strand = "original"
    
    if k_mer_rc < min_k_mer:
      min_k_mer = k_mer_rc
      min_pos = j
      min_strand = "reverse_complement"
  
  return (min_k_mer, i + min_pos + offset, min_strand, k, w)

def find_minimizers(sequence, k, w):
  n = len(sequence)
  indices = range(n - w - k + 2)
  args = [(sequence, k, w, i,0) for i in indices]

  with mp.Pool() as pool:
    minimizers = pool.map(find_window_minimizer, args)

  return minimizers



# if __name__ == "__main__":
#     sequence = "ACGTTGCAACGTTGCA"
#     k = 3
#     w = 5
#     minimizers = find_minimizers(sequence, k, w)
#     for minimizer in minimizers:
#       print(f"Minimizer: {minimizer[0]}, Position: {minimizer[1]}, Origin: {minimizer[2]}, k: {minimizer[3]}, w: {minimizer[4]}")
   