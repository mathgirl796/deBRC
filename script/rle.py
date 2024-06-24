import numpy as np
import sys
from collections import Counter, OrderedDict
import huffman
from tqdm import tqdm

def err_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

def rle(x):
    """Find runs of consecutive items in an array."""

    # ensure array
    x = np.asanyarray(x)
    if x.ndim != 1:
        raise ValueError('only 1D array supported')
    n = x.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(x[:-1], x[1:], out=loc_run_start[1:])
        run_starts = np.nonzero(loc_run_start)[0]

        # find run values
        run_values = x[loc_run_start]

        # find run lengths
        run_lengths = np.diff(np.append(run_starts, n))

        return run_values, run_lengths, run_starts
    
def dna_compress(values, lengths):
    dna_format_lengths = []
    flat_dna_format_lengths = []
    split_bits_num = 2
    mod_factor = 2 ** split_bits_num
    dna_format_lengths = None
    err_print('txt to bits to chars')
    while np.any(lengths != 0):
        remain = lengths % mod_factor
        remain[lengths <= 0] = -1
        # err_print(remain)
        if dna_format_lengths is None:
            dna_format_lengths = remain.reshape((1, -1))
        else:
            dna_format_lengths = np.concatenate((dna_format_lengths, remain.reshape((1, -1))), axis=0)
        lengths = lengths // 4
    dna_format_lengths[(dna_format_lengths >= 0)] += 256
    # err_print(dna_format_lengths.shape, dna_format_lengths)
    err_print('building huffman tree')
    source_code = np.concatenate((dna_format_lengths.T, values.reshape((-1, 1))), axis=1) # 每行代表lengths中的一个数，由低位到高位，后面都是-1，表示没用，最后连接上value
    recipe = list(zip(*np.unique(source_code[(source_code != -1)], return_counts=True)))
    codebook = huffman.codebook(recipe)
    codebook[-1] = ''
    err_print(codebook)
    # err_print('huffman encoding')
    out = [codebook[x] for x in tqdm(source_code.flatten(), desc='huffman_encoding') if x != -1]
    err_print('concatenating')
    out = ''.join(out)
    # err_print(out)
    return out, codebook
    
    
if __name__ == '__main__':
    s = '#FFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFF::FFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,:FFFFFFFFFFFFFFFFFFFFFFF,FF:FFF:F:FFF:FFFFFFFFF'
    ar = np.frombuffer(s.encode('ascii'), dtype=np.uint8)
    value, length, _ = rle(np.array(list(ar)))
    dna_compress(value, length)
