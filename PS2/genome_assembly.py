import numpy as np
import argparse
from pprint import pprint
from collections import deque, defaultdict
import random
from copy import deepcopy

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
    parser = argparse.ArgumentParser(description="CS482 PS2 Question 4 solution a29asgha")
    parser.add_argument('--input', type=str, required=True, help="Input FASTA file containing DNA sequences")
    parser.add_argument('--output', type=str, default="assembled.fna", help="Output path for the result of the program")
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

def get_substrings(s, k):
    for i in range(len(s) - k + 1):
        yield s[i:i+k]

def reverse_complement(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def build_graph(seqs: list[str], k= 0) -> list[tuple[str, str]]:
    res = []
    for seq in seqs:
        for s in get_substrings(seq, len(seq) - k):
            res.append((s[:-1], s[1:]))
    
    return res
    
def count_components_bfs(adj_list):
    visited = set()
    components = 0

    def bfs(start):
        """Iterative BFS to visit all nodes in a component."""
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for neighbor in adj_list.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

    for node in adj_list:
        if node not in visited:
            visited.add(node)
            bfs(node)
            components += 1

    return components


def remove_until_one_tip_and_sink(adj_list):
    # helper function used inside of the main loop
    def find_incoming_edges(adj_list: dict):
    
        incoming_edges = defaultdict(int)
        for node, neighbors in adj_list.items():
            for neighbor in neighbors:
                incoming_edges[neighbor] += 1
        return incoming_edges
    tips = None
    sinks = None
    while True:
        incoming_edges = find_incoming_edges(adj_list)
        outgoing_edges = {node: len(neighbors) for node, neighbors in adj_list.items()}
        tips = {node for node in adj_list if incoming_edges[node] == 0}
        sinks = {node for node in adj_list if outgoing_edges[node] == 0}

        # If we have exactly one tip and one sink, stop
        if len(tips) == 1 and len(sinks) == 1:
            break

        if len(tips) > 1:
            for tip in list(tips)[:-1]:
                del adj_list[tip]

        if len(sinks) > 1:
            for sink in list(sinks)[:-1]:
                for node in adj_list:
                    adj_list[node] = [n for n in adj_list[node] if n != sink]

    return adj_list, tips.pop(), sinks.pop()


def correct_errors(seqs: list, k):
    k_mers = {}
    for h, v in seqs:
        for mer in get_substrings(v, k):
            if mer not in k_mers:
                k_mers[mer] = 0    
            k_mers[mer] += 1
    
    lower_freqs = [key for key, value in k_mers.items() if value <= 5]
    higher_freqs = [key for key, value in k_mers.items() if value > 5]
    
    def hamming_distance(seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    
    replacements = {}
    for u in lower_freqs:
        best = len(v)
        for v in higher_freqs:
            d = hamming_distance(u, v)
            if (d < best):
                replacements[u] = v
                best = d

    def correct_sequence(seq):
        """Replace error k-mers in a sequence with corrected versions."""
        corrected_seq = list(seq)
        
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if kmer in replacements:
                corrected_kmer = replacements[kmer]
                corrected_seq[i : i + k] = corrected_kmer  # Replace part of sequence

        return "".join(corrected_seq)
    
    corrected = [correct_sequence(seq) for h, seq in seqs]
    return corrected


def eulerian_cycle(graph, start_node):
    def find_cycle(start_node):
        cycle = []
        node = start_node
        while graph[node]:
            next_node = random.choice(graph[node])
            cycle.append((node, next_node))
            graph[node].remove(next_node)
            node = next_node
        return cycle
    graph_copy = {node: edges[:] for node, edges in graph.items()}
    # start_node = next((node for node, edges in graph_copy.items() if edges), None)
    
    if start_node is None:
        return [] 
    
    cycle = find_cycle(start_node)
    return cycle

def main():
    parser = make_parser()
    args = parser.parse_args()
    seqs = read_fasta_file(args.input)
    print("starting error correction")
    corrected = correct_errors(seqs, 21)
    
    print("initial graph constructed, stats:")
    edges = build_graph(corrected, 80)
    
    nodes = {}
    for u, v in edges:
        if u not in nodes:
            nodes[u] = []
        if v not in nodes:
            nodes[v] = []
        nodes[u].append(v)
    print("number of nodes: ", len(nodes))
    print("number of edges: ", sum([len(v) for k, v in nodes.items()]))
    print("number of components: ", count_components_bfs(nodes))
    print("Cleaning the graph by removing sinks and tips:")
    nodes, t, s = remove_until_one_tip_and_sink(nodes)
    nodes[s].append(t)
    print("Graph cleaning done; new stats:")
    print("number of nodes: ", len(nodes))
    print("number of edges: ", sum([len(v) for k, v in nodes.items()]))
    print("number of components: ", count_components_bfs(nodes))
    
    cycle = eulerian_cycle(deepcopy(nodes), t)
    print(len(cycle))
    
    idx = 0
    for i in range(len(cycle)):
        if (cycle[i][0] == t):
            idx = i
    cycle = cycle[idx:] + cycle[:idx]
    retVal = ""
    for v, u in cycle:
        retVal += v[0]
    retVal += cycle[-1][0][1:]
    with open(args.output, mode="w") as f:
        f.write(">assembled \n")
        f.write(retVal)
    
    
    print("Program done!")


if __name__ == "__main__":
    main()