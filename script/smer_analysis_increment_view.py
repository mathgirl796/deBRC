from matplotlib import pyplot as plt
import json
import numpy as np
import os
import argparse
import sys

def err_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputFilePath', type=str, required=True)
args = parser.parse_args()

filepath = args.inputFilePath

err_print('start load file')
with open(filepath, 'r') as f:
    di = json.load(f)
countBaseVar=di['countBaseVar'] if 'countBaseVar' in di else False

err_print('start compute title')
x = np.array(di['x'], dtype=np.uint64)
y = np.array(di['y'], dtype=np.uint64)
average = np.round(np.dot(x, y.T) / np.sum(y), 2)
cumulative_delta = np.dot(x, y.T)
title = f"file={os.path.basename(di['filepath'])}, k={di['k']}, delta_type_num={len(x)}, countBaseVar={countBaseVar}\n"\
        +f"kmer_num={int(np.sum(y))+1}, cumulative_delta={cumulative_delta}\n"\
        +f"average_delta={average}, min_delta={np.min(x)}, max_delta={np.max(x)}"
err_print(title)
plt.figure(figsize=(12,10))
plt.title(title)
plt.xscale('linear' if countBaseVar else 'log')
plt.xlabel('delta')
plt.yscale('log')
plt.ylabel('occurance')

err_print('start draw figure')
if countBaseVar:
    plt.bar(x, y)
else:
    sample_num = np.min((cumulative_delta, 10000))
    sample_index = np.random.choice(np.arange(len(x)), sample_num, p=y/np.sum(y))
    err_print(f"sample {len(sample_index)} data")
    title += f"sample_num={sample_num}"
    plt.title(title)
    x = x[sample_index]
    y = y[sample_index]
    plt.bar(x, y)

err_print('start save figure')
plt.savefig(filepath + ".png")