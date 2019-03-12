#!/usr/bin/env python

import argparse
import os

argparser = argparse.ArgumentParser(description = 'Merges SumHer output files.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_list_file', required = True, help = 'Input file with a list of SumHer output files (one per line).')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')


def read_sumher_output(in_file):
    pheno1, pheno2, suffix = os.path.split(in_file)[-1].split('.')
    pheno1, pheno2 = map(lambda x: x.replace('__', '.'), [pheno1, pheno2])
    with open(in_file, 'rt') as ifile:
        header = ifile.readline()
        if not header: return {}
        header = header.rstrip().split()
        if not all(h in header for h in ['Component', 'Value', 'SD']): return {}
        t_header = ['Pheno1', 'Pheno2']
        t_row = [pheno1, pheno2]
        for line in ifile:
            line = line.rstrip()
            if not line: continue
            columns = line.split()
            if len(columns) != len(header): continue
            columns = dict(zip(header, columns))
            t_header.append(columns['Component'])
            t_row.append(columns['Value'])
            t_header.append(columns['Component'] + '_SD')
            t_row.append(columns['SD'])
        return { 'header': t_header, 'record': t_row }


if __name__ == '__main__':
    args = argparser.parse_args()
    header = False
    in_files = []
    with open(args.in_list_file, 'rt') as ifile:
        for line in ifile:
            line = line.strip()
            if not line: continue
            in_files.append(line)
    with open(args.out_file, 'wt') as ofile:
        for in_file in in_files:
            corr = read_sumher_output(in_file)
            if not corr: continue
            if not header:
                ofile.write('{}\n'.format('\t'.join(corr['header'])))
                header = True
            ofile.write('{}\n'.format('\t'.join(corr['record'])))
