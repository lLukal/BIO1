def find_window_minimizer(args):
  sequence, k, w, i,strand = args
  window = sequence[i:i+w+k-1]
  min_k_mer = window[:k]
  min_pos = 0
  
  for j in range(1, w):
    k_mer = window[j:j+k]
    
    if k_mer < min_k_mer:
      min_k_mer = k_mer
      min_pos = j
    
  return (min_k_mer, i + min_pos, strand, k, w)

def find_minimizers(sequence, k, w,strand):
  n = len(sequence)
  indices = range(n - w - k + 2)
  windows = [(sequence, k, w, i,strand) for i in indices]

  minimizers = []
  for window in windows:
    minimizers.append(find_window_minimizer(window))

  return minimizers
   