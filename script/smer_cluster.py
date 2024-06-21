import struct
from collections import Counter, defaultdict
import json
import numpy as np
import sys
import argparse


def err_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

filepath = "/home/user/duanran/repo/deBRC/deBRC/experiment/single_genomes/data/chm13.draft_v1.1.f1_assembly_v2_genbank/k21/chm13.draft_v1.1.f1_assembly_v2_genbank.smer"
cluster_len_list = [10, 9, 8, 7, 6]

input_file = open(filepath, "rb")
kmerLen, totalKmerNum = struct.unpack('<IQ', input_file.read(12))
err_print(f"kmerLen\tkmerNum\n{kmerLen}\t{totalKmerNum}")
input_file.close()

# batch_kmer_num = 128500000
batch_kmer_num = np.iinfo(np.uint32).max
progress = 0
total_progress = int(np.ceil((totalKmerNum - 1) / batch_kmer_num))
result = defaultdict(list)
err_print(f"progress:{progress}/{total_progress}")
while progress * batch_kmer_num < totalKmerNum:
    a = np.fromfile(filepath, dtype=np.uint64, count=batch_kmer_num, offset=12 + batch_kmer_num * 8 * progress)
    a = np.right_shift(a, (kmerLen-10)*2)
    for cluster_len in cluster_len_list:
        err_print(f"computing cluster_len {cluster_len}")
        unique_elements, counts = np.unique(a, return_counts=True)
        # err_print(unique_elements, counts)
        result[cluster_len].append(counts)
        a = np.right_shift(a, 2)
    progress += 1
    err_print(f"progress:{progress}/{total_progress}")

for cluster_len in cluster_len_list:
    result[cluster_len] = np.concatenate(result[cluster_len])
    err_print(f"saving cluster {cluster_len}, shape: {result[cluster_len].shape}")
    np.save(filepath+f".{cluster_len}.npy", result[cluster_len])
