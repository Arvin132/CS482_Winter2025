import sys
import argparse
from enum import Enum
from fasta import read_fasta_file  # Custom module to read FASTA files

# Function to define and parse command-line arguments
def make_parser():
    parser = argparse.ArgumentParser(description="CS482 PS1 Question 3 solution a29asgha")
    parser.add_argument('--db', type=str, required=True, help="FASTA file containing the sequences to search for")
    parser.add_argument('--query', type=str, required=True, help="FASTA file containing a single sequence acting as the query")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    parser.add_argument('--k', type=int, default=11, help="Length of each k-mer")
    return parser

# Enum class representing possible paths in the alignment process
class PATH(Enum):
    DIAG_MATCH = 0
    RIGHT = 1
    DOWN = 2
    DIAG_MISMATCH = 3
    BEGIN = 4

# Function to build a dictionary of k-mers and their positions in the sequence
def build_k_mers(seq: str, k: int = 5) -> dict[str, list[int]]:
    retVal = {}
    i = 0
    while i < len(seq) - k:
        mer = seq[i: i + k]
        if retVal.get(mer) is None:  # Initialize list if k-mer not already in the dictionary
            retVal[mer] = []
        retVal[mer].append(i)  # Store the starting index of the k-mer
        i += 1
    return retVal

# Function to calculate a similarity score between two sequences
def calculate_score(s1: str, s2: str, match_score: float = 1.0, mismatch_penalty: float = -1.0, gap_penalty: float = -1.0) -> float:
    if len(s2) > len(s1):  # Ensure s1 is the longer sequence
        calculate_score(s2, s1)
    
    score = 0.0
    for i in range(len(s1)):
        if i >= len(s2):  # Handle gaps when sequences are of different lengths
            score -= gap_penalty
        elif s1[i] == s2[i]:  # Match
            score += match_score
        elif s1[i] == '-' or s2[i] == '-':  # Gap
            score -= gap_penalty
        else:  # Mismatch
            score -= mismatch_penalty
    return score

# Function to find the best match for a query sequence in a database sequence
def search_best_match(seq1: str, seq2: str, mers: dict[str, list[int]]) -> tuple[float, int, int, int]:
    hsp = 0.0  # Highest scoring pair
    match_loc = None
    k = None
    left_extension = 0
    right_extension = 0

    for mer, lst in mers.items():
        if k is None:
            k = len(mer)  # Set k to the length of the k-mers
        i = 0
        while i < len(seq1) - k:
            score = calculate_score(mer, seq1[i: i + k])
            
            # Variables to track the best left and right extension
            best_right_end = None
            best_left_end = None

            if score > k - 1:  # Check if initial score is significant
                best_addon_left = 0.0
                best_addon_right = 0.0
                
                # Extend matches to the left
                for loc in lst:
                    right = i + k
                    left = i - 1
                    query_left = loc - 1
                    query_right = loc + k
                    addon = 0.0
                    misses = 0
                    
                    # Extend to the left until mismatch threshold is reached
                    while left >= 0 and query_left >= 0 and misses <= 2:
                        addon += 1 if seq1[left] == seq2[query_left] else -1
                        misses += 1 if seq1[left] != seq2[query_left] else 0
                        if addon > best_addon_left:
                            best_addon_left = addon
                            best_left_end = left
                        left -= 1
                        query_left -= 1

                    # Extend to the right until mismatch threshold is reached
                    addon = 0.0
                    misses = 0
                    while right < len(seq1) and query_right < len(seq2) and misses <= 2:
                        addon += 1 if seq1[right] == seq2[query_right] else -1
                        misses += 1 if seq1[right] != seq2[query_right] else 0
                        if addon > best_addon_right:
                            best_addon_right = addon
                            best_right_end = right
                        right += 1
                        query_right += 1
                
                # Update total score
                score += best_addon_left + best_addon_right
            
            # Update highest scoring pair if a better score is found
            if score > hsp:
                if best_left_end is not None:
                    left_extension = i - best_left_end
                if best_right_end is not None:
                    right_extension = best_right_end - i + k - 1
                hsp = score
                match_loc = i
            i += 1

    return hsp, match_loc, left_extension, right_extension

# Main function to parse arguments, process data, and find the best match
def main():
    parser = make_parser()
    args = parser.parse_args()
    db = read_fasta_file(args.db)  # Read database sequences
    query: str = read_fasta_file(args.query)[0][1]  # Read query sequence
    k = args.k
    mers = build_k_mers(query, k)  # Build k-mer dictionary

    best_score = 0.0
    best_header = ""
    best_idx = None
    best_loc = None
    best_le = 0
    best_re = 0

    # Iterate through database entries to find the best match
    for idx, (head, content) in enumerate(db):
        if args.verbose: 
            print("Current entry being checked:", head)
        score, loc, le, re = search_best_match(content, query, mers)
        if score > best_score:
            best_score = score
            best_header = head
            best_loc = loc
            best_idx = idx
            best_le = le
            best_re = re

    # Print the best match and its details
    print(f"Best Match in the database: {best_header}")
    print(f"Match score: {best_score}")
    print("Matched sequence:")
    print(db[best_idx][1][best_loc - best_le: best_loc + best_re + 1])

if __name__ == "__main__":
    main()
