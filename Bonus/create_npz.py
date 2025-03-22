import numpy as np
from itertools import product

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
seq = "AAAATAGCTANNNNTAGACGTNNNGTAGTC"
k = 3
res = find_frequency(seq, k)
f = np.array([v for k, v in res.items()])
f = f / f.sum()
name = "11_multiple_comp"
np.savez(name, f)
print(len(seq) - seq.count("N"), k)
with open(name + "_output.txt", mode="w+") as f:
    f.write(seq)
    f.write(f"\n L={len(seq) - seq.count("N")}, k={k}")