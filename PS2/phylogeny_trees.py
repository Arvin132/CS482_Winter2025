from Bio import Entrez, SeqIO
import argparse
from itertools import product
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo import draw_ascii

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
    parser = argparse.ArgumentParser(description="CS482 PS2 Question 3 solution a29asgha")
    parser.add_argument('--input', type=str, required=True, help="Input file containing accession ids of organisms")
    parser.add_argument('--email', type=str, default="asgharianarvin@gmail.com", help="email used to connect to NCBI server")
    parser.add_argument('--test', type=str, default="candidate.fas", help="the test candidate")
    parser.add_argument('--k', default=4, type=int, help="size of each mer")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

def get_substrings(s, k):
    for i in range(len(s) - k + 1):
        yield s[i:i+k]

def read_ids(filename: str) -> list[str]:
    retVal = []
    
    with open(filename, mode="r+") as f:
        line = f.readline()
        while(line.strip() != "" and line is not None):
            line = line[line.find(":") + 1:]
            tokens = line.split(',')
            for token in tokens:
                retVal.append(token.strip())
            line = f.readline()
    
    return retVal

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

def download_genome(accession_id, email):
    try:
        # Fetch genome data from NCBI
        Entrez.email = email
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = handle.read()
        return record[record.find('\n') + 1:].strip()
    except Exception as e:
        print("Error: ", e)  

def create_distance_matrix(vecs, distance_func) -> list[list[float]]:
    D = np.zeros(shape=(len(vecs), len(vecs)))
    retVal = []
    
    
    for i in range(len(vecs)):
        cur = []
        for j in range(i + 1):
            if (i != j):
                cur.append(distance_func(vecs[i], vecs[j]))
            else:
                cur.append(0)
        retVal.append(cur)
        
    return retVal

     

def main():
    parser = make_parser()
    args = parser.parse_args()
    ids = read_ids(args.input)
    seqs = []
    print("Ids read, now Downloading sequences")
    for id in ids:
        seqs.append(download_genome(id, args.email))
    
    vecs = []
    print("Sequences reading done! now calculating k-mer frequencies")
    for seq in seqs:
        result = find_frequency(seq, args.k)
        vecs.append(list(result.values()))
    vecs = np.array(vecs)
    vecs = vecs / np.linalg.norm(vecs, ord=1) # normalize with L1 norm
    
    
    cosine_similarity = lambda v1, v2: float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    pearson_corr = lambda v1, v2: float(np.cov(v1, v2, bias=True)[0, 1] / (np.std(v1) * np.std(v2)))
    
    methods = [
        ("Euclidean",  lambda i, j: float(np.linalg.norm(i - j))),
        ("Cosine",  lambda v1, v2: float(1 - cosine_similarity(v1, v2))),
        ("Pearson",  lambda v1, v2: float(1 - pearson_corr(v1, v2)))
    ]
    
    
    for name, d_func in methods:
        print(f"Now creating upgma tree for {name} distance!")
        dist = create_distance_matrix(vecs, d_func)
        distance_matrix = DistanceMatrix(ids, dist)
        constructor = DistanceTreeConstructor()
        upgma_tree = constructor.upgma(distance_matrix)
        draw_ascii(upgma_tree, open(name + ".nwk", mode="w+"))
    
    print("Now adding candiadate and rerunning the analysis!")
    seq = read_fasta_file(args.test)[0][1]
    ids.append("test")
    vecs = list(vecs)
    vecs.append(list(find_frequency(seq, args.k).values()))
    vecs = np.array(vecs)
    
    for name, d_func in methods:
        print(f"Now creating upgma tree for {name} distance!")
        dist = create_distance_matrix(vecs, d_func)
        distance_matrix = DistanceMatrix(ids, dist)
        constructor = DistanceTreeConstructor()
        upgma_tree = constructor.upgma(distance_matrix)
        draw_ascii(upgma_tree, open(name + "_wcandidate.nwk", mode="w+"))
    
    

if __name__ == "__main__":
    main()