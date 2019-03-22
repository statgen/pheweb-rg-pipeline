import argparse
import json
import os

argparser = argparse.ArgumentParser(description = 'Chunks pipline into independent runs based on JSON phenotypes file. Note: creates multiple directories named "chunk_*".')
argparser.add_argument('-i', '--in-json', metavar = 'file', dest = 'in_json', required = True, help = 'JSON phenotypes file.')
argparser.add_argument('-s', '--size', metavar = 'number', dest = 'chunk_size', type = int, required = True, help = 'Max number of phenotypes in a chunk.')

if __name__ == '__main__':
   args = argparser.parse_args()
   if args.chunk_size <= 1:
       print('Max number of phenotypes in a chunk must be greater than 1.')
       exit(1)
   with open(args.in_json) as f_in:
      all_data = json.load(f_in)
   all_pheno_list = []
   compute_pheno_list = []
   for c, i in enumerate(range(0, len(all_data), args.chunk_size), 1):
      dirname = f'chunk_{c}'
      os.mkdir(dirname)
      compute_pheno_list = all_data[i:i + args.chunk_size]
      all_pheno_list.extend(compute_pheno_list)
      with open(os.path.join(dirname, 'pheno-list.json'), 'w') as f_out:
         json.dump(all_pheno_list, f_out, indent = 4)
      with open(os.path.join(dirname, 'compute-pheno-list.json'), 'w') as f_out:
         json.dump(compute_pheno_list, f_out, indent = 4)
