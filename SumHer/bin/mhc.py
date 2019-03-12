#!/usr/bin/env python

import argparse
import gzip

argparser = argparse.ArgumentParser(description = 'Lists MHC variants in *.bim files.')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'in_files', nargs = '+', required = True, help = 'Input *.bim file(s).')
argparser.add_argument('-b', '--genome-build', metavar = 'name', dest = 'build', choices = ['GRCh37', 'GRCh38'], required = True, help = 'Human genome build version: GRCh37, GRCh38.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')

mhc = {
        'GRCh37': { 'chrom': '6', 'start': 28477797, 'stop': 33448354  } ,
        'GRCh38': { 'chrom': '6', 'start': 28510120, 'stop': 33480577  }
}

if __name__ == '__main__':
    args = argparser.parse_args()
    mhc_chrom = mhc[args.build]['chrom']
    mhc_start = mhc[args.build]['start']
    mhc_stop = mhc[args.build]['stop']
    with open(args.out_file, 'wt') as ofile:
        for in_file in args.in_files:
            with open(in_file, 'rt') as ifile:
                for line in ifile:
                    chrom, name, morgans, bp, a1, a2  = line.rstrip().split()
                    if chrom == mhc_chrom:
                        bp = int(bp)
                        if bp >= mhc_start and bp <= mhc_stop:
                            ofile.write('{}\n'.format(name))
