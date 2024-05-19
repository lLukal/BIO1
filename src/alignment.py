import numpy as np

# Global Alignment - Needleman-Wunsch
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    # Initialize the scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m+1, n+1))
    
    # Initialize the first row and column
    for i in range(1, m+1):
        score_matrix[i][0] = i * gap
    for j in range(1, n+1):
        score_matrix[0][j] = j * gap
    
    # Fill the scoring matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score_matrix[i-1][j] + gap
            insert = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match_score, delete, insert)
    
    # The alignment score is the value in the bottom-right cell
    alignment_score = score_matrix[m][n]
    return score_matrix, alignment_score


# Local Alignment - Smith-Waterman
def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    # Initialize the scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m+1, n+1))
    
    # Fill the scoring matrix
    max_score = 0
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score_matrix[i-1][j] + gap
            insert = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match_score, delete, insert, 0)
            max_score = max(max_score, score_matrix[i][j])
    
    # The alignment score is the maximum value in the matrix
    alignment_score = max_score
    return score_matrix, alignment_score


# Semi-Global Alignment
def semi_global(seq1, seq2, match=1, mismatch=-1, gap=-1):
    # Initialize the scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m+1, n+1))
    
    # Initialize the first row and column with gaps
    for i in range(1, m+1):
        score_matrix[i][0] = 0
    for j in range(1, n+1):
        score_matrix[0][j] = 0
    
    # Fill the scoring matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score_matrix[i-1][j] + gap
            insert = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match_score, delete, insert)
    
    # Allow free end gaps
    max_score = max(max(score_matrix[m]), max(score_matrix[:,n]))
    
    # The alignment score is the maximum value in the last row or column
    alignment_score = max_score
    return score_matrix, alignment_score


# # Example sequences
# seq1 = "GATTACA"
# seq2 = "GCATGCU"

# # Compute the alignments
# nw_matrix, nw_score = needleman_wunsch(seq1, seq2)
# sw_matrix, sw_score = smith_waterman(seq1, seq2)
# sg_matrix, sg_score = semi_global(seq1, seq2)

# # Print the scoring matrices and the alignment scores
# print("Needleman-Wunsch (Global Alignment):")
# print(nw_matrix)
# print(f"Alignment score: {nw_score}")

# print("\nSmith-Waterman (Local Alignment):")
# print(sw_matrix)
# print(f"Alignment score: {sw_score}")

# print("\nSemi-Global Alignment:")
# print(sg_matrix)
# print(f"Alignment score: {sg_score}")