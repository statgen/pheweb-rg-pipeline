#!/usr/bin/env python

import argparse
import os

argparser = argparse.ArgumentParser(description = 'Intersects files preserving the order. Note: suitable only for very short lines (e.g. rsId).')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_files', nargs = '+', required = True, help = 'Input files.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')


if __name__ == '__main__':
    args = argparser.parse_args()
    sizes = []
    for in_file in args.in_files:
        sizes.append((os.stat(in_file).st_size, in_file))
    sizes.sort(key = lambda x: x[0])
    names = dict()
    with open(sizes[0][1], 'rt') as ifile:
        for line in ifile:
            line = line.rstrip()
            if line: names[line] = 1
    for size, in_file in sizes[1:]:
        with open(in_file, 'rt') as ifile:
            for line in ifile:
                line = line.rstrip()
                if line in names: names[line] += 1
    with open(args.out_file, 'wt') as ofile:
        with open(sizes[0][1], 'rt') as ifile:
            for line in ifile:
                line = line.rstrip()
                if names[line] == len(sizes):
                    ofile.write('{}\n'.format(line))
