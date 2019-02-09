import argparse
import json

argparser = argparse.ArgumentParser(description = 'Generates pheno-list.json file from tab-delimited file. Warning: will overwrite existing pheno-list.json file.')
argparser.add_argument('-i', '--in-tsv', metavar = 'file', dest = 'in_tsv', required = True, help = 'Input tab-delimited file. No header. Column order: full path to summary stat file, phenocode, number of cases, number of controls.')

if __name__ == '__main__':
   args = argparser.parse_args()
   data = []
   with open(args.in_tsv, 'r') as f_in:
      for line in f_in:
         columns = line.rstrip().split('\t')
         if any(len(c.strip()) == 0 for c in columns):
             continue
         if len(columns) < 4:
             continue
         path = columns[0]
         phenocode = columns[1]
         cases = columns[2]
         controls = columns[3]
         try:
             cases_int = int(cases)
         except:
             raise Exception('Non integer value \'{}\' found in cases column.'.format(cases))
         try:
             controls_int = int(controls)
         except:
             raise Exception('Non integer value \'{}\' found in controls column.'.format(controls))
         data.append({
            'assoc_files': [path],
            'phenocode': phenocode, 
            'num_cases': cases_int, 
            'num_controls': controls_int
         })
   with open('pheno-list.json', 'w') as f_out:
      json.dump(data, f_out, indent = 4)
