import matplotlib.pyplot as plt

def plot_mapped_genome(results, reference_length):
    fig, ax = plt.subplots(figsize=(15, 5))
    
    # Plot the reference genome as a horizontal line
    ax.plot([0, reference_length], [0, 0], color='blue', label='Reference Genome')
    
    # Plot each fragment as a horizontal line positioned below the reference genome
    for idx, result in enumerate(results):
        fragment_result, ref_name, q_begin, q_end, t_begin, t_end, aligned_seq1, aligned_seq2, alignment_score, cigar = result
        y_position = -(idx + 1)
        ax.plot([t_begin, t_end], [y_position, y_position], marker='o', label=f'Fragment {idx+1}')

    # Adding labels and title
    ax.set_xlabel('Reference Genome Position')
    ax.set_yticks([0] + [-(i+1) for i in range(len(results))])
    ax.set_yticklabels(['Reference'] + [f'Fragment {i+1}' for i in range(len(results))])
    ax.set_title('Mapped Genome Visualization')
    ax.legend(loc='upper right')
    ax.grid(True)

    plt.savefig("output.png")
    # Display the plot
    plt.show()

# # Example results array
# results = [
#     ("fragment1", "reference_genome", 0, 50, 100, 150, "AGCTAGCTAG", "AGCTAGCTAG", 60),
#     ("fragment2", "reference_genome", 0, 30, 300, 330, "GCTAGCTAGC", "GCTAGCTAGC", 50),
#     ("fragment3", "reference_genome", 0, 20, 500, 520, "TTCGATCGAT", "TTCGATCGAT", 45)
# ]
# reference_length = 1000

# # Plot the mapped genome
# plot_mapped_genome(results, reference_length)
