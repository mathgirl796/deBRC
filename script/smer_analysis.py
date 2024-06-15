import struct
from collections import Counter, defaultdict
import json
import numpy as np
import sys
import argparse

def err_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--smerFilePath', type=str, required=True)
parser.add_argument('-b', '--countBaseVar', action='store_true')
args = parser.parse_args()

filepath = args.smerFilePath
countBaseVar = args.countBaseVar
if len(sys.argv) > 1 and sys.argv[1] == '1':
    countBaseVar = True


input_file = open(filepath, "rb")
kmerLen, totalKmerNum = struct.unpack('<IQ', input_file.read(12))
err_print(f"kmerLen\tkmerNum\n{kmerLen}\t{totalKmerNum}")
input_file.close()

first_kmer = np.fromfile(filepath, dtype=np.uint64, count=1, offset=12)[0]
last_kmer = np.fromfile(filepath, dtype=np.uint64, count=1, offset=12 + 8*(totalKmerNum - 1))[0]
average_distance = (last_kmer - first_kmer) / (totalKmerNum - 1)
err_print(f"first_kmer:{bin(first_kmer)}, last_kmer:{bin(last_kmer)}, average_distance:{average_distance}")

batch_kmer_num = 128500000
progress = 0
total_progress = int(np.ceil((totalKmerNum - 1) / batch_kmer_num))
result = defaultdict(int)
err_print(f"progress:{progress}/{total_progress}")
while progress * batch_kmer_num < totalKmerNum:
    a = np.fromfile(filepath, dtype=np.uint64, count=batch_kmer_num + 1, offset=12 + batch_kmer_num * 8 * progress)
    deltas = a[1:] - a[:-1]
    if countBaseVar:
        diff_nums = np.zeros_like(deltas)
        diff_tail_base_num = 1
        while not np.all(deltas == 0):
            deltas = np.right_shift(deltas, 2)
            # err_print(np.count_nonzero(deltas))
            diff_nums[(deltas == 0) & (diff_nums == 0)] = diff_tail_base_num
            diff_tail_base_num += 1
        deltas = diff_nums
    ct = Counter(deltas.tolist())
    for key, value in ct.items():
        result[key] += value
    progress += 1
    err_print(f"progress:{progress}/{total_progress}")

x, y = zip(*result.items())
output = {'k':kmerLen, 'filepath':filepath, 'countBaseVar':countBaseVar, 'x': x, 'y': y}
# err_print(output)
with open(filepath + ".baseNumVar.json" if countBaseVar else filepath + ".kmerDistance.json", "w") as f:
    json.dump(output, f)
err_print("done dump result to json file")

