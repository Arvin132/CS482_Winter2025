# Importing necessary modules for command-line arguments and enum definitions
import sys
import argparse
from enum import Enum
from fasta import read_fasta_file  # Custom function to read sequences from a FASTA file

# Creates an argument parser for CLI input
def make_parser():
    parser = argparse.ArgumentParser(description="CS482 PS1 Question 1 solution a29asgha")
    parser.add_argument('--input', type=str, default="", help="Input FASTA file")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

# Enum to represent possible backtracking paths in the alignment matrix
class PATH(Enum):
    DIAG_MATCH = 0
    RIGHT = 1
    DOWN = 2
    DIAG_MISMATCH = 3

# Computes the alignment score (S) and backtracking matrix (B)
def find_best_alignment(seq1: str, seq2: str, match_score=1.0, gap_penalty=-1.0, mismatch_penalty=0.0):
    n = len(seq1)  # Length of sequence 1
    m = len(seq2)  # Length of sequence 2

    # Initialize scoring and backtracking matrices
    S = [[0.0 for i in range(m + 1)] for i in range(n + 1)]
    B = [[None for i in range(m + 1)] for i in range(n + 1)]

    # Fill the first column and row with gap penalties
    for i in range(n + 1):
        S[i][0] = i * gap_penalty
        B[i][0] = PATH.DOWN  # Down indicates a gap in seq2
    for j in range(m + 1):
        S[0][j] = j * gap_penalty
        B[0][j] = PATH.RIGHT  # Right indicates a gap in seq1

    B[0][0] = PATH.DIAG_MATCH  # Initialization for the starting cell

    # Fill the scoring and backtracking matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            add = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            match = S[i - 1][j - 1] + add
            delete = S[i - 1][j] + gap_penalty
            insert = S[i][j - 1] + gap_penalty

            # Determine the best score and path
            ma = max(match, delete, insert)
            if ma == match:
                B[i][j] = PATH.DIAG_MATCH if add == match_score else PATH.DIAG_MISMATCH
            elif ma == delete:
                B[i][j] = PATH.DOWN
            else:
                B[i][j] = PATH.RIGHT
            S[i][j] = ma  # Update the score matrix

    return S, B  # Return scoring and backtracking matrices

# Reconstructs the best alignment from the backtracking matrix
def create_best_alignment(seq1: str, seq2: str, B):
    i = len(seq1)
    j = len(seq2)
    retVal1 = ''  # Aligned sequence 1
    retVal2 = ''  # Aligned sequence 2

    while True:
        if i == 0 and j == 0:  # Reached the start of the matrix
            return retVal1, retVal2
        if B[i][j] in (PATH.DIAG_MATCH, PATH.DIAG_MISMATCH):
            retVal1 = seq1[i - 1] + retVal1
            retVal2 = seq2[j - 1] + retVal2
            i -= 1
            j -= 1
        elif B[i][j] == PATH.DOWN:
            retVal1 = seq1[i - 1] + retVal1
            retVal2 = '-' + retVal2
            i -= 1
        else:  # PATH.RIGHT
            retVal1 = '-' + retVal1
            retVal2 = seq2[j - 1] + retVal2
            j -= 1

# Calculates the edit distance using the backtracking matrix
def calculate_edit_distance(B):
    i = len(B) - 1
    j = len(B[0]) - 1
    retVal = 0  # Edit distance

    while True:
        if i == 0 and j == 0:  # Reached the start
            return retVal
        if B[i][j] in (PATH.DIAG_MATCH, PATH.DIAG_MISMATCH):
            if B[i][j] == PATH.DIAG_MISMATCH:
                retVal += 1  # Mismatch counts as an edit
            i -= 1
            j -= 1
        elif B[i][j] == PATH.DOWN:
            retVal += 1  # Gap in seq2
            i -= 1
        else:  # PATH.RIGHT
            retVal += 1  # Gap in seq1
            j -= 1

# Main function for parsing arguments and executing alignment
def main():
    parser = make_parser()
    args = parser.parse_args()
    fasta = read_fasta_file(args.input)  # Read sequences from FASTA
    seq1 = fasta[0][1]
    seq2 = fasta[1][1]

    S, B = find_best_alignment(seq1, seq2)

    if args.verbose:  # Display detailed output if verbose mode is enabled
        import numpy as np
        print(np.array(S))  # Print scoring matrix
        print(np.array(B))  # Print backtracking matrix
        s1, s2 = create_best_alignment(seq1, seq2, B)
        print("Best Alignments")
        print("Sequence 1: ", s1)
        print("Sequence 2: ", s2)

    print("Distance of 2 sequences: ")
    print(calculate_edit_distance(B))  # Print the edit distance

# Entry point
if __name__ == "__main__":
    main()
