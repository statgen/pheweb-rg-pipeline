#!/usr/bin/env python

import argparse
import gzip

argparser = argparse.ArgumentParser(description = 'Formats summary statistics')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'in_file', required = True, help = 'Input file with summary statistics, compressed using gzip/bgzip.')
argparser.add_argument('-o', '--out', metavar = 'prefix', dest = 'out_prefix', required = False, help = 'Prefix for output files: <prefix>.stats, <prefix>.nonamb, <prefix>.big.')
argparser.add_argument('-s', '--snp', metavar = 'column', dest = 'snp', default = 'ID', help = 'Column name for SNP.')
argparser.add_argument('-ea', '--effect-allele', metavar = 'column', dest = 'ea', default = 'ALT', help = 'Column name for effect allele.')
argparser.add_argument('-oa', '--other-allele', metavar = 'column', dest = 'oa', default = 'REF', help = 'Column name for other allele.')
argparser.add_argument('-e', '--effect-size', metavar = 'column', dest = 'beta', default = 'beta', help = 'Column name for effect size (log odds).')
argparser.add_argument('-se', '--se-effect', metavar = 'column', dest = 'se', default = 'sebeta', help = 'Column name for standard error of effect size.')
argparser.add_argument('-n', '--num-samples', metavar = 'column', nargs = '+', dest = 'n', default = ['num_cases', 'num_controls'], help = 'Column name(s) for number of samples. If multiple columns are specified, they are summped up to get total sample size.')


if __name__ == '__main__':
    args = argparser.parse_args()
    required_columns = [ args.snp, args.ea, args.oa, args.beta, args.se ] + args.n
    stats_filename = f'{args.out_prefix}.stats'
    pred_filename = f'{args.out_prefix}.nonamb'
    le_pred_filename = f'{args.out_prefix}.big'
    unique_names = set()
    with gzip.open(args.in_file, 'rt') as ifile, \
            open(stats_filename, 'w') as stats_ofile, \
            open(pred_filename, 'w') as pred_ofile, \
            open(le_pred_filename, 'w') as le_pred_ofile:
        iheader = ifile.readline()
        if not iheader:
            raise Exception(f'Empty header in {args.in_file}.')
        iheader = iheader.rstrip().split()
        for c in required_columns:
            if c not in iheader:
                raise Exception(f'Missing {c} column in {args.in_file}.')
        stats_ofile.write('Predictor A1 A2 Direction Stat n\n')
        for i, iline in enumerate(ifile, 2):
            iline = iline.rstrip().split()
            if len(iline) != len(iheader):
                raise Exception(f'Number of columns on line {i} does not match header.')
            record = dict(zip(iheader, iline))
            if record[args.snp] in unique_names: # skip if ID is duplicated
                continue
            unique_names.add(record[args.snp])
            a1 = record[args.ea].upper()
            a2 = record[args.oa].upper()
            try:
                n = sum(map(int, (record[c] for c in args.n)))
                stat = (float(record[args.beta]) / float(record[args.se])) ** 2
            except: # catch possible type conversion or division by 0 errors
                continue
            if stat > n / 99.0:
                le_pred_ofile.write('{}\n'.format(record[args.snp]))
            if not (a1 in {'A', 'T'} and a2 in {'A', 'T'}) and not (a1 in {'C', 'G'} and a2 in {'C', 'G'}):
                pred_ofile.write('{}\n'.format(record[args.snp]))
            stats_ofile.write('{} {} {} {} {} {}\n'.format(record[args.snp], a1, a2, record[args.beta], stat, n))
