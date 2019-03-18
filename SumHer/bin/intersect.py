#!/usr/bin/env python

import argparse
from collections import Counter

argparser = argparse.ArgumentParser(description = 'Intersects files preserving the order. Note: suitable only for very short lines (e.g. rsId).')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_files', nargs = '+', required = True, help = 'Input files.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')


if __name__ == '__main__':
    args = argparser.parse_args()
    n_files = len(args.in_files)
    names = Counter()
    for in_file in args.in_files:
        with open(in_file, 'rb') as ifile:
            names.update(ifile)
    with open(args.out_file, 'wb') as ofile:
        with open(args.in_files[0], 'rb') as ifile:
            for line in ifile:
                if names[line] == n_files:
                    ofile.write(line)
