import sys
from itertools import product
import numpy as np
from time import time
import argparse


def err_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

def bin2fa(args):
    of = open(args.output, 'w')
    motif_map = {i:''.join(x) for i,x in enumerate(product('ACGT', repeat=4))}
    for fn in args.inputs:
        t1 = time()
        with open(fn, 'rb') as f:
            data = f.read(args.fragment_size)
            fragment_id = 0
            while len(data) > 0:
                im = np.frombuffer(data, dtype=np.uint8)
                err_print(fn, type(data), len(data), im.dtype, im.shape)
                motifs = np.vectorize(motif_map.get)(im.reshape(-1))
                of.write(f'>{fragment_id}|{fn}\n{''.join(motifs.tolist())}\n')
                err_print(f'time: {time()-t1}s', motifs.dtype)
                data = f.read(args.fragment_size)
                fragment_id += 1
    of.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, required=True, help='output file path')
    parser.add_argument('-f', '--fragment_size', type=int, required=False, default=None, help='cut each file into fragments')
    parser.add_argument("inputs", nargs="*", help="input file list")
    args = parser.parse_args()
    bin2fa(args)