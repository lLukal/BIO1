
def reverse_complement(sequence):
    """
    Function to find the reverse complement of a DNA sequence.

    Parameters:
    sequence (str): The DNA sequence to find the reverse complement of.

    Returns:
    str: The reverse complement of the input DNA sequence.

    """
    # Dictionary to store the complement of each nucleotide
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    # Find the reverse complement by mapping each nucleotide to its complement, reversing the sequence, and joining the nucleotides
    return ''.join(complement[base] for base in reversed(sequence))

def kmer_to_uint(kmer):
    """
    Function to convert a kmer to an unsigned integer.

    Parameters:
    kmer (str): The kmer to convert.

    Returns:
    int: The unsigned integer representation of the kmer.
    """
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    binary = ''.join(mapping[nucleotide] for nucleotide in kmer)
    return int(binary, 2)
     


def Minimize(sequence,sequence_len,kmer_len, window_len):
    """
    Function to find minimizers in a given sequence.

    Parameters:
    sequence (str): The DNA sequence to find minimizers in.
    kmer_len (int): The length of the k-mers to consider.
    window_len (int): The length of the window to consider for finding minimizers.

    Returns:
    list: A list of tuples, where each tuple contains the hash of a minimizer, its position in the sequence, and a boolean indicating whether it comes from the original strand (True) or the reverse complement (False).

    """

    # List to store the minimizers
    minimizers = []

    for i in range(0, sequence_len - window_len + 1):

        # Extract the current window from the sequence
        window = sequence[i:i+window_len]

        # Find the lexicographically smallest kmer in the current window
        min_kmer = min(window[j:j+kmer_len] for j in range(window_len - kmer_len + 1))

        # Compute the reverse complement of the smallest kmer
        min_kmer_reverse = str(reverse_complement(min_kmer))

        if min_kmer < min_kmer_reverse:
            minimizers.append((kmer_to_uint(min_kmer),i,True))
        else:
            minimizers.append((kmer_to_uint(min_kmer_reverse),i,False))

    # Return the list of minimizers
    return minimizers


def main():
    sequence = "ATGCGATCGTACGTA"
    sequence_len = len(sequence)
    kmer_len = 3
    window_len = 5

    minimizers = Minimize(sequence, sequence_len, kmer_len, window_len)

    for minimizer, pos, is_original_strand in minimizers:
        strand = "original" if is_original_strand else "reverse complement"
        print(f"Minimizer: {minimizer}, Position: {pos}, Strand: {strand}")

if __name__ == "__main__":
        main()