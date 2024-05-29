import numpy as np

# Global Alignment - Needleman-Wunsch
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
  m, n = len(seq1), len(seq2)
  score_matrix = np.zeros((m+1, n+1))
  traceback_matrix = np.zeros((m+1, n+1), dtype=np.int8)
  
  # Initialize the scoring matrix and traceback matrix
  for i in range(1, m+1):
    score_matrix[i][0] = i * gap
    traceback_matrix[i][0] = 1
  for j in range(1, n+1):
    score_matrix[0][j] = j * gap
    traceback_matrix[0][j] = 2
  
  for i in range(1, m+1):
    for j in range(1, n+1):
      match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
      delete = score_matrix[i-1][j] + gap
      insert = score_matrix[i][j-1] + gap
      score_matrix[i][j], traceback_matrix[i][j] = max((match_score, 0), (delete, 1), (insert, 2))
  
  # Traceback
  aligned_seq1, aligned_seq2 = '', ''
  i, j = m, n
  while i > 0 or j > 0:
    if traceback_matrix[i][j] == 0:
      aligned_seq1 = seq1[i-1] + aligned_seq1
      aligned_seq2 = seq2[j-1] + aligned_seq2
      i -= 1
      j -= 1
    elif traceback_matrix[i][j] == 1:
      aligned_seq1 = seq1[i-1] + aligned_seq1
      aligned_seq2 = '-' + aligned_seq2
      i -= 1
    else:
      aligned_seq1 = '-' + aligned_seq1
      aligned_seq2 = seq2[j-1] + aligned_seq2
      j -= 1
  
  alignment_score = score_matrix[m][n]
  return score_matrix, aligned_seq1, aligned_seq2, alignment_score, ''


# Local Alignment - Smith-Waterman
def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
  m, n = len(seq1), len(seq2)
  score_matrix = np.zeros((m+1, n+1))
  traceback_matrix = np.zeros((m+1, n+1), dtype=np.int8)
  
  max_i, max_j = 0, 0
  max_score = 0
  
  for i in range(1, m+1):
    for j in range(1, n+1):
      match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
      delete = score_matrix[i-1][j] + gap
      insert = score_matrix[i][j-1] + gap
      score_matrix[i][j], traceback_matrix[i][j] = max((0, -1), (match_score, 0), (delete, 1), (insert, 2))
      if score_matrix[i][j] >= max_score:
          max_score = score_matrix[i][j]
          max_i, max_j = i, j
  
  # Traceback
  aligned_seq1, aligned_seq2 = '', ''
  i, j = max_i, max_j
  while i > 0 and j > 0 and score_matrix[i][j] > 0:
    if traceback_matrix[i][j] == 0:
      aligned_seq1 = seq1[i-1] + aligned_seq1
      aligned_seq2 = seq2[j-1] + aligned_seq2
      i -= 1
      j -= 1
    elif traceback_matrix[i][j] == 1:
      aligned_seq1 = seq1[i-1] + aligned_seq1
      aligned_seq2 = '-' + aligned_seq2
      i -= 1
    else:
      aligned_seq1 = '-' + aligned_seq1
      aligned_seq2 = seq2[j-1] + aligned_seq2
      j -= 1
  
  alignment_score = max_score
  return score_matrix, aligned_seq1, aligned_seq2, alignment_score, ''


# Semi-Global Alignment
def semi_global(seq1, seq2, match=1, mismatch=-1, gap=-1):
  m, n = len(seq1), len(seq2)
  score_matrix = np.zeros((m+1, n+1))
  traceback_matrix = np.zeros((m+1, n+1), dtype=np.int8)
  
  for i in range(1, m+1):
    score_matrix[i][0] = 0
    traceback_matrix[i][0] = 1
  for j in range(1, n+1):
    score_matrix[0][j] = 0
    traceback_matrix[0][j] = 2
  
  for i in range(1, m+1):
    for j in range(1, n+1):
      match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
      delete = score_matrix[i-1][j] + gap
      insert = score_matrix[i][j-1] + gap
      score_matrix[i][j], traceback_matrix[i][j] = max((match_score, 0), (delete, 1), (insert, 2))
  
  # Traceback
  aligned_seq1, aligned_seq2 = '', ''
  
  i, j = m, np.argmax(score_matrix[m]) if np.max(score_matrix[m]) > np.max(score_matrix[:, n]) else np.argmax(score_matrix[:, n])
  while i > 0 or j > 0:
    if traceback_matrix[i][j] == 0:
      aligned_seq1 = seq1[i-1] + aligned_seq1
      aligned_seq2 = seq2[j-1] + aligned_seq2
      i -= 1
      j -= 1
    elif traceback_matrix[i][j] == 1:
      aligned_seq1 = seq1[i-1] + aligned_seq1
      aligned_seq2 = '-' + aligned_seq2
      i -= 1
    else:
      aligned_seq1 = '-' + aligned_seq1
      aligned_seq2 = seq2[j-1] + aligned_seq2
      j -= 1
  
  alignment_score = max(np.max(score_matrix[m]), np.max(score_matrix[:, n]))
  return score_matrix, aligned_seq1, aligned_seq2, alignment_score, ''


# # Example sequences
# seq1 = "GATTACA"
# seq2 = "GCATGCU"

# # Compute the alignments and get the scores
# nw_matrix, nw_aligned_seq1, nw_aligned_seq2, nw_score = needleman_wunsch(seq1, seq2)
# sw_matrix, sw_aligned_seq1, sw_aligned_seq2, sw_score = smith_waterman(seq1, seq2)
# sg_matrix, sg_aligned_seq1, sg_aligned_seq2, sg_score = semi_global(seq1, seq2)

# # Print the scoring matrices and the alignment scores
# print("Needleman-Wunsch (Global Alignment):")
# print(nw_matrix)
# print(f"Aligned Sequences:\n{nw_aligned_seq1}\n{nw_aligned_seq2}")
# print(f"Alignment score: {nw_score}")

# print("\nSmith-Waterman (Local Alignment):")
# print(sw_matrix)
# print(f"Aligned Sequences:\n{sw_aligned_seq1}\n{sw_aligned_seq2}")
# print(f"Alignment score: {sw_score}")

# print("\nSemi-Global Alignment:")
# print(sg_matrix)
# print(f"Aligned Sequences:\n{sg_aligned_seq1}\n{sg_aligned_seq2}")
# print(f"Alignment score: {sg_score}")