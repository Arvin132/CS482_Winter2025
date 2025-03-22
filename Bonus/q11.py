import numpy as np
import argparse
from copy import deepcopy
import random

def make_parser():
    parser = argparse.ArgumentParser(description="CS482 Midterm bonus solution a29asgha")
    parser.add_argument('--L', type=int, required=True, help="Length of the wanted string to be constructed")
    parser.add_argument('--k', type=int, required=True, help="Length of each k-mer")
    parser.add_argument('--p', type=str, required=True, help="Path to the probability vector")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser


class Graph:
    def __init__(self):
        self.nodes = {}
    
    def __getitem__(self, node):
        return self.nodes[node]
    
    def add_node(self, u):
        if (u not in self.nodes):
            self.nodes[u] = []
            
    def add_edge(self, u, v):
        self.add_node(u)
        self.add_node(v)
        self.nodes[u].append(v)
    
    def pop_edge(self, u ,v):
        self.nodes[u].remove(v)
    
    def __repr__(self):
        return self.nodes.__repr__()  
    
    def in_degree(self, u):
        retVal = 0
        for node in self.nodes:
            retVal += self.nodes[node].count(u)
        return retVal
    
    def __len__(self):
        return len(self.nodes)
    
    def out_degree(self, u):
        return len(self.nodes[u])
    

def build_debruijn_graph(mers: list[str]) -> Graph:
            
    retVal = Graph()
    for mer in mers:
        prefix = mer[:-1]
        suffix = mer[1:]
        retVal.add_edge(prefix, suffix)
    return retVal
     
def find_eulerian_cycle(graph: Graph, start) -> list:
    curpath = [start]
    retVal = []
    while(len(curpath) > 0):
        curNode = curpath[-1] 
        if (len(graph[curNode]) > 0):
            nextNode = graph[curNode][-1]
            graph[curNode].pop()
            curpath.append(nextNode)
        else:
            retVal.append(curpath[-1])
            curpath.pop()
    
    retVal = retVal[::-1]
    return retVal

def build_kmer_array(freq: np.ndarray, k) -> list[str]:
    def get_index_to_kmer(idx: int, k) -> str:
        retVal = ""
        trans = {0 : "A", 1: "C", 2: "G", 3: "T"}
        for _ in range(k):
            retVal += trans[idx % 4]
            idx = idx //4
        return retVal[::-1]
    
    mers = []
    for i, v in enumerate(freq):
        for _ in range(v):
            mers.append(get_index_to_kmer(i, k))
    
    return mers
    
def main():
    parser = make_parser()
    args = parser.parse_args()
   
    k = args.k
    L = args.L
    p = next(iter((np.load(args.p)).values()))
    p_unorm = np.rint(p * (L - k + 1)).astype(int)
    mers = build_kmer_array(p_unorm, k) # return list[str]
    if (k == 1):
        random
        print("".join(random.sample(mers, len(mers))))
        return
    graph = build_debruijn_graph(mers)
    
    # Find the start and end of the strongly connected graph
    start = None
    end = None
    for node in graph.nodes:
        if graph.in_degree(node) != graph.out_degree(node):
            if graph.in_degree(node) < graph.out_degree(node):
                start = node
            else:
                end = node
    
    circular_seq = False
    if (start == None and end == None):
        start = end = next(iter(graph.nodes))
        circular_seq = True
    
    aug_graph = deepcopy(graph)
    if not circular_seq: aug_graph.add_edge(end, start)
    cycle = find_eulerian_cycle(aug_graph, start)
    cycle.pop()
    
    created_seq = ""
    for node in cycle:
        created_seq += node[0]
    created_seq += cycle[-1][1:]
    if (circular_seq):
        created_seq += cycle[0][-1]
    
    if (args.verbose):
        print("Given graph is:")
        print(graph)
        print("Found Cycle:")
        print(cycle)
    
    if (len(created_seq) > L):
        created_seq = created_seq[:L]
    if (len(created_seq) < L):
        created_seq = created_seq + (''.join(random.choice(['A','G','T','C']) for _ in range(L - len(created_seq))))
    
    print(created_seq, end="")
    
    
    


if __name__ == "__main__":
    main()