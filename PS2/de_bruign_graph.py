import numpy as np
import argparse

# Creates an argument parser for CLI input
def make_parser():
    parser = argparse.ArgumentParser(description="CS482 PS2 Question 2 solution a29asgha")
    parser.add_argument('--input', type=str, required=True, help="Input file containing DNA sequences")
    parser.add_argument('--output', type=str, default="deBruijn.txt", help="Output path for the result of the program")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

def read_seq(input_file: str) -> list[str]:
    retVal = []
    with open(input_file, mode='r') as f:
        
        seq = f.readline().strip()
        while(seq is not None and seq != ""):
            retVal.append(seq)
            seq = f.readline().strip()  
    return retVal

def reverse_complement(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def get_substrings(s, k):
    for i in range(len(s) - k + 1):
        yield s[i:i+k]

def build_graph(seqs: list[str]) -> set[tuple[str, str]]:
    res = set()
    for seq in seqs:
        rc = reverse_complement(seq)
        res.add((seq[:-1], seq[1:]))
        res.add((rc[:-1], rc[1:]))
    
    return res
    
    
        
        

def main():
    parser = make_parser()
    args = parser.parse_args()
    seqs = read_seq(args.input)
    res = build_graph(seqs)
    with open(args.output, mode='w+') as f:
        for g in sorted(res):
            f.write(str(g) + '\n')
    
    

if __name__ == "__main__":
    main()