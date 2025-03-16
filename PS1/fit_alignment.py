import sys
import argparse
from enum import Enum
from fasta import read_fasta_file

def make_parser():
    # Creates a command-line argument parser for input and verbosity options
    parser = argparse.ArgumentParser(description="CS482 PS1 Question 2 solution a29asgha")
    parser.add_argument('--input', type=str, default="", required=True, help="Input FASTA file")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

class PATH(Enum):
    # Enum for possible traceback directions in alignment
    DIAG_MATCH = 0  # Diagonal match
    RIGHT = 1       # Right (gap in seq1)
    DOWN = 2        # Down (gap in seq2)
    DIAG_MISMATCH = 3  # Diagonal mismatch
    BEGIN = 4       # Start point for traceback

def fit_best_alignment(seq1: str, seq2: str, match_score=1.0, gap_penalty=-1.0, mismatch_penalty=-1.0) -> tuple[list[list[int]], list[list[PATH]]]:
    # Performs local alignment on two sequences
    n = len(seq1)
    m = len(seq2)
    
    # Initialize scoring and traceback matrices
    S = [[0.0 for i in range(m + 1)] for i in range(n + 1)]
    B = [[None for i in range(m + 1)] for i in range(n + 1)]
    
    # Fill first column and row with gap penalties
    for i in range(n + 1):
        S[i][0] = i * gap_penalty
        B[i][0] = PATH.DOWN
    for j in range(m + 1):
        S[0][j] = j * gap_penalty
        B[0][j] = PATH.RIGHT
    B[0][0] = PATH.BEGIN
    
    # Fill the rest of the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            add = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            match = S[i - 1][j - 1] + add
            delete = S[i - 1][j] + gap_penalty
            insert = S[i][j - 1] + gap_penalty
            
            ma = max(match, delete, insert, 0)  # Use 0 for local alignment

            # Determine traceback direction
            if ma == 0:
                B[i][j] = PATH.BEGIN
            elif ma == match:
                if add == match_score:
                    B[i][j] = PATH.DIAG_MATCH
                else:
                    B[i][j] = PATH.DIAG_MISMATCH
            elif ma == delete:
                B[i][j] = PATH.DOWN
            else:
                B[i][j] = PATH.RIGHT
            S[i][j] = ma

    return S, B

def create_fit(seq1: str, seq2: str, S, B) -> tuple[str, str]:
    # Builds the best alignment by tracing back from the highest scoring cell
    i = 0
    j = 0
    retVal1 = ''
    retVal2 = ''
    
    # Find the cell with the highest score
    for x in range(len(seq1)):
        for y in range(len(seq2)):
            if S[x][y] >= S[i][j]:
                i = x
                j = y

    # Traceback to build the alignment
    while True:
        if B[i][j] == PATH.BEGIN:
            return retVal1, retVal2
        if B[i][j] == PATH.DIAG_MATCH or B[i][j] == PATH.DIAG_MISMATCH:
            retVal1 = seq1[i - 1] + retVal1
            retVal2 = seq2[j - 1] + retVal2
            i -= 1
            j -= 1
        elif B[i][j] == PATH.DOWN:
            retVal1 = seq1[i - 1] + retVal1
            retVal2 = '-' + retVal2
            i -= 1
        else:  # B[i][j] == PATH.RIGHT
            retVal1 = '-' + retVal1
            retVal2 = seq2[j - 1] + retVal2
            j -= 1

def main():
    # Main function to parse input, compute alignments, and display results
    parser = make_parser()
    args = parser.parse_args()
    fasta = read_fasta_file(args.input)
    seq1 = fasta[0][1]
    seq2 = fasta[1][1]
    
    # Compute scoring and traceback matrices
    S, B = fit_best_alignment(seq1, seq2)
    
    # Display matrices if verbose flag is set
    if args.verbose:
        import numpy as np
        print(np.array(S))
        print(np.array(B))
    
    # Create and display the best alignment
    s1, s2 = create_fit(seq1, seq2, S, B)
    print("Best Alignments")
    print("Sequence 1: ", s1)
    print("Sequence 2: ", s2)
    
    # Find and display the best alignment score
    i, j = 0, 0
    for x in range(len(seq1)):
        for y in range(len(seq2)):
            if S[x][y] >= S[i][j]:
                i = x
                j = y
    print(f"\nBest Score: {S[i][j]}")

if __name__ == "__main__":
    main()
