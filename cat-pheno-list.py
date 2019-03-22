import argparse
import json

argparser = argparse.ArgumentParser(description = 'Concatenates JSON files with phenotypes.')
argparser.add_argument('-i', '--in-json', metavar = 'file', dest = 'in_jsons', nargs = '+', required = True, help = 'List of JSON phenotype files.')
argparser.add_argument('-o', '--out-json', metavar = 'file', dest = 'out_json', required = True, help = 'Output JSON phenotype file.')

if __name__ == '__main__':
   args = argparser.parse_args()
   data = []
   for in_json in args.in_jsons:
      with open(in_json) as f_in:
         data.extend(json.load(f_in))
   with open(args.out_json, 'w') as f_out:
      json.dump(data, f_out, indent = 4)
