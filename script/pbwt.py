import sys, os
import numpy as np
from collections import Counter
import cv2
import argparse
from tqdm import tqdm

def trans(args):
    filepath, outputpath = args.input, args.output
    start, end = args.start, args.end
    lines = []
    with open(filepath) as f:
        for i in range(start):
            f.readline()
        for i in range(end - start):
            line = f.readline()
            if len(line) > 0:
                lines.append(line)
            else:
                break
    quals = [x.strip() for x in lines]
    qual_lens = [len(x) for x in quals]
    qual_distribution = Counter(qual_lens)
    print(f'quan_len distribution: {qual_distribution}')

    bmp = np.zeros((len(quals), np.max(list(qual_distribution.keys()))), dtype=np.uint8)
    for i, qual in enumerate(tqdm(quals, desc='trans')):
        b = np.frombuffer(qual.encode('ascii'), dtype=np.uint8)
        b.astype(np.uint8)
        bmp[i] = b

    if outputpath is None:
        outputpath = filepath + '.bmp'
    cv2.imwrite(outputpath, bmp)

def pbwt_encode_core(bmp):
    # encode pbwt
    pbwt = np.zeros_like(bmp)
    pbwt[:,0] = bmp[:,0]
    order = np.arange(bmp.shape[0])
    for j in tqdm(list(range(1, bmp.shape[1])), desc='encode'):
        index = np.argsort(bmp[:,j-1][order], kind='stable') # retrive and sort new column
        order = order[index]
        pbwt[:,j] = bmp[:,j][order] # save result
    return pbwt

def pbwt_encode(args):
    filepath, outputpath, transpose = args.input, args.output, args.transpose
    bmp = cv2.imread(filepath, cv2.IMREAD_GRAYSCALE)
    pbwt = pbwt_encode_core(bmp)
    if outputpath is None:
        if transpose:
            outputpath = filepath + '.pbwt.T.webp'
        else:
            outputpath = filepath + '.pbwt.webp'
    cv2.imwrite(outputpath, pbwt.T if transpose else pbwt)

def pbwt_decode_core(pbwt):
    # decode pbwt
    rbmp = np.zeros_like(pbwt)
    rbmp[:,0] = pbwt[:,0]
    order = np.arange(pbwt.shape[0])
    for j in tqdm(list(range(1, pbwt.shape[1])), desc='decode'):
        index = np.argsort(rbmp[:,j-1][order], kind='stable')
        order = order[index]
        rbmp[:,j][order] = pbwt[:,j]
    return rbmp

def pbwt_decode(args):
    filepath, outputpath, transpose = args.input, args.output, args.transpose
    pbwt = cv2.imread(filepath, cv2.IMREAD_GRAYSCALE)
    if transpose:
        pbwt = pbwt.T
    rbmp = pbwt_decode_core(pbwt)
    if outputpath is None:
        outputpath = filepath + '.bmp'
    cv2.imwrite(outputpath, rbmp)

def to_webp(args):
    filepath, outputpath = args.input, args.output
    img = cv2.imread(filepath)
    if outputpath is None:
        outputpath = filepath + '.webp'
    cv2.imwrite(outputpath, img)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='pbwt')
    subparsers = parser.add_subparsers(required=True)
    parser_trans = subparsers.add_parser('trans', help='transform a txt to bmp')
    parser_trans.add_argument('-i', '--input', type=str, required=True, help='input file path')
    parser_trans.add_argument('-o', '--output', type=str, required=False, default=None, help='output file path')
    parser_trans.add_argument('-t', '--transpose', action='store_true', help='transpose img before output')
    parser_trans.add_argument('-s', '--start', type=int, required=False, default=0, help='start line num')
    parser_trans.add_argument('-e', '--end', type=int, required=False, default=214748647, help='end line num')
    parser_encode = subparsers.add_parser('encode', help='encode grayscale bmp to pbwt.T.webp format')
    parser_encode.add_argument('-i', '--input', type=str, required=True, help='input file path')
    parser_encode.add_argument('-o', '--output', type=str, required=False, default=None, help='output file path')
    parser_encode.add_argument('-t', '--transpose', action='store_true', help='transpose pbwt before output')
    parser_encode.add_argument('-w', '--webp', action='store_true', help='output webp format (default bmp)')
    parser_decode = subparsers.add_parser('decode', help='decode a pbwt.T.webp format file to bmp')
    parser_decode.add_argument('-i', '--input', type=str, required=True, help='input file path')
    parser_decode.add_argument('-o', '--output', type=str, required=False, default=None, help='output file path')
    parser_decode.add_argument('-t', '--transpose', action='store_true', help='transpose pbwt before output')
    parser_webp = subparsers.add_parser('webp', help='convert an image to webp format')
    parser_webp.add_argument('-i', '--input', type=str, required=True, help='input file path')
    parser_webp.add_argument('-o', '--output', type=str, required=False, default=None, help='output file path')

    parser_trans.set_defaults(func=trans)
    parser_encode.set_defaults(func=pbwt_encode)
    parser_decode.set_defaults(func=pbwt_decode)
    parser_webp.set_defaults(func=to_webp)

    args = parser.parse_args()
    args.func(args)