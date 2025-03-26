import argparse
import collections


def pattern_matching_with_suffix_array(Q, suffixes ):
    
    left, right = 0, len(suffixes)
    # First Binary Search for L_Q (First occurrence)
    while left < right:
        mid = (left + right) // 2
        if suffixes[mid][0] < Q: # Move right if suffix is smaller
            left = mid + 1
        else:
            right = mid # Move left otherwise
    first = left
    right = len(suffixes)
    # Second Binary Search for U_Q (Last occurrence)
    while left < right:
        mid = (left + right) // 2
        if suffixes[mid][0].startswith(Q): # If Q is a prefix, move right
            left = mid + 1
        else: # Otherwise, move left
            right = mid
    last = right-1 # Adjust to return last valid match
    
    if first > last:
        return None
    else:
        return (first, last)

def BWT(given: str):
    rotations = []
    for i in range(len(given)):
        rotations.append(given[i:] + given[:i])
    print(rotations)
    rotations = sorted(rotations)
    print(rotations)
    return "".join([v[-1] for v in rotations])

def BWTi(given: str):
    acc = []
    for ch in given:
        acc.append(ch)
    while(len(acc[0]) != len(given)):
        acc = sorted(acc)
        for i in range(len(acc)):
            acc[i] = given[i] + acc[i]
        
    return next(x for x in acc if x.endswith('$'))

def make_parser():
    parser = argparse.ArgumentParser(description="CS482 PS3 Q1 solution a29asgha")
    parser.add_argument("input", type=str, help="Path to the input file")
    parser.add_argument("--output", default="indices.txt", type=str, help="Path to the output file")
    parser.add_argument('--verbose', '-v', action='store_true', help="Enable verbose output")
    return parser

def main():
    parser = make_parser()
    args = parser.parse_args()
    seqs = []
    with open(args.input, mode="r+") as f:
        r = f.readline().strip()
        while(len(r) != 0 and r is not None):
            seqs.append(r)
            r = f.readline().strip()
    
    ref: str = seqs[0]
    queries = seqs[1:]
    
    # ref = BWTi(ref)
    suffix_array = sorted([(ref[i:], i) for i in range(len(ref))], key=lambda kv: kv[0])
    matches = []
    for query in queries:
        res = pattern_matching_with_suffix_array(query, suffix_array)
        if (res != None):
            for idx in range(res[0], res[1] + 1):
                matches.append(suffix_array[idx][1])
    # matches = sorted(matches)
    
    with open(args.output, mode="w+") as f:
        for match in matches:
            f.write(str(match))
            if (args.verbose): print(match)
            f.write("\n")
    
    


if __name__ == "__main__":
    main()