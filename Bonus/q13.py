import numpy as np
import argparse
from collections import defaultdict

def reverse_complement(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def read_fasta_file(file: str) -> list[tuple[str, str]]:
    
    """
    Reads a FASTA file and extracts the headers and sequences.

    Parameters:
    -----------
    file : str
        Path to the FASTA file to be read.

    Returns:
    --------
    list[tuple[str, str]]
        Each tuple contains:
        - header
        - The corresponding sequence
    
    Raises:
    -------
    ValueError
        If the file is not in a valid FASTA format (e.g., missing '>' at the start).
    FileNotFoundError
        If the specified file does not exist.

    Notes:
    ------
    - Assumes the FASTA file does not contain any multi-line headers.
    - The sequence is treated as a single string without newlines.
    """
    
    retVal = []
    
    with open(file, mode='r+') as f:
        ch = ch = f.read(1)
        while(True): 
            if(ch == ""): break
            if (ch != '>'): raise ValueError("Bad Input FASTA file")
            
            header = ""
            ch = f.read(1)
            while(ch != "\n"):
                header += ch
                ch = f.read(1)
                
            content = ""
            ch = f.read(1)
            while(ch != ">" and ch != ""):
                if (ch != '\n'):
                    content += ch
                ch = f.read(1)
            
            retVal.append((header, content))
    
    return retVal


def make_parser():
    parser = argparse.ArgumentParser(description="CS482 PS2 Question 2 solution a29asgha")
    parser.add_argument('--reads', type=str, default="test_0.fasta", help="Input file containings reads of the program")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

def main():
    parser = make_parser()
    args = parser.parse_args()
    reads = read_fasta_file(args.reads)
    genomes = [v for h, v in reads]
    originals = set(genomes)
    table = defaultdict(int)
    for genome in genomes:
        table[genome] += 1
        table[reverse_complement(genome)] += 1
    
    errors = [k for k, v in table.items() if v <= 1 and k in originals]
    replacement = []
    
    for error in errors:
        for i in range(len(error)):
            is_break = False
            for rep in ["A", "C", "G", "T"]:
                if (error[i] == rep):
                    continue
                possible = error[0:i] + rep + error[i + 1:]
                if (possible in table and table[possible] > 1):
                    replacement.append(possible)
                    is_break = True
                    break
            if (is_break):
                break
        
    if (len(replacement) != len(errors)):
        raise ValueError("For some sequences, no matching replacement was found, raising error")
    
    for err, replace in zip(errors, replacement):
        print(err, "->", replace, sep="")


if __name__ == "__main__" :
    main()