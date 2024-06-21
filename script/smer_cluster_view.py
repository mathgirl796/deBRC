import struct
from collections import Counter, defaultdict
import json
import numpy as np
import sys
import argparse
from matplotlib import pyplot as plt

def err_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

filepath = "/home/user/duanran/repo/deBRC/deBRC/experiment/single_genomes/data/chm13.draft_v1.1.f1_assembly_v2_genbank/k21/chm13.draft_v1.1.f1_assembly_v2_genbank.smer"
cluster_len_list = [10, 9, 8, 7, 6]
quantiles=[0.2, 0.8]

input_file = open(filepath, "rb")
kmerLen, totalKmerNum = struct.unpack('<IQ', input_file.read(12))
err_print(f"kmerLen\tkmerNum\n{kmerLen}\t{totalKmerNum}")
input_file.close()

datas = []
for cluster_len in cluster_len_list:
    a = np.load(filepath+f".{cluster_len}.npy")
    err_print(a.shape)
    datas.append(a)
err_print('done read data')


err_print('drawing violin')
plt.figure(figsize=(12,10))
plt.title(f"quantiles={quantiles}")
plt.xscale('linear')
plt.xlabel('l')
plt.yscale('log')
plt.ylabel('cluster_size')
plt.xticks(range(len(cluster_len_list) + 1), labels=[0]+cluster_len_list)
plt.grid(True, axis='y', which='major', linestyle='--', color='gray', alpha=0.5)
plt.violinplot(datas, showmedians=True, quantiles=[quantiles for a in range(len(cluster_len_list))])
plt.savefig("/home/user/duanran/repo/deBRC/deBRC/shit/violin.png")

err_print('drawing plot')
plt.show()
plt.figure(figsize=(30,10))
plt.xscale('log')
plt.xlabel('cluster_size')
plt.yscale('log')
plt.ylabel('cluster_num')
for i in range(len(cluster_len_list)):
    unique_items, counts = np.unique(datas[i], return_counts=True)
    plt.plot(unique_items, counts, linewidth=2, label=cluster_len_list[i], alpha=0.5)
plt.grid(True, axis='y', which='major', linestyle='--', color='gray', alpha=0.5)
plt.legend()
plt.savefig("/home/user/duanran/repo/deBRC/deBRC/shit/plot.png")