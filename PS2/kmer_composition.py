import argparse
from itertools import product


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



# Creates an argument parser for CLI input
def make_parser():
    parser = argparse.ArgumentParser(description="CS482 PS2 Question 1 solution a29asgha")
    parser.add_argument('--k', type=int, help="size of each mer", required=True)
    parser.add_argument('--input', type=str, required=True, help="Input Fasta file containing a DNA sequence")
    parser.add_argument('--output', type=str, default="k-mer.txt", help="Output path for the result of the program")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

def get_substrings(s, k):
    for i in range(len(s) - k + 1):
        yield s[i:i+k]

def find_frequency(seq: str, k: int) -> dict[str, int]:
    mers = {}
    newSeq = []
    valid_chars = {'A', 'G', 'C', 'T', 'N'}
    
    # filter the sequence to get rid of bad characters
    for ch in seq:
        if (ch in valid_chars):
            newSeq.append(ch)
        else:
            newSeq.append('N')
    newSeq = "".join(newSeq)
    # initilize the k-mer array
    for mer in product(['A', 'C', 'G', 'T'], repeat=k):
        val = "".join(mer)
        mers[val] = 0
    
    for sub in get_substrings(newSeq, k):
        if (sub in mers):
            mers[sub] += 1
    
    return mers
    
    
    

def main():
    parser = make_parser()
    args = parser.parse_args()
    seq = read_fasta_file(args.input)[0][1]
    res = find_frequency(seq, args.k)
    # print the output to the needed file
    with open(args.output, mode='w+') as f:
        for c in res.values():
            f.write(str(c) + '\n')
    
    

if __name__ == "__main__":
    main()