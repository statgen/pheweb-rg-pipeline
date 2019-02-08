import json
import os
import itertools


configfile: "config.json"


def get_all_phenocodes(json_file):
   with open(json_file, 'r') as f_in:
      traits = json.load(f_in)
      for trait in traits:
         yield trait['phenocode']   


def read_corr(filename):
   with open(filename, 'r') as f_in:
      while True:
         line = f_in.readline()
         if not line: break
         if line.startswith("Summary of Genetic Correlation Results"): break
      header = f_in.readline()
      if not header: return {}
      header = header.rstrip().split()
      record = f_in.readline()
      if not record: return {}
      record = dict(zip(header, record.rstrip().split()))
      for p in ['p1', 'p2']:
         record[p] = re.sub(r'\.sumstats\.gz$', '', os.path.split(record[p])[-1])
      return {'header': header, 'record': record}


def write_merged(data, filename):
   with open(filename, 'w') as f_out:
      header = data[0]['header']
      f_out.write('{}\n'.format('\t'.join(header)))
      for d in data:
         f_out.write('{}\n'.format('\t'.join([d['record'][h] for h in header])))


rule all:
   input: "result/ALL.RG.txt"


rule merge:
   input: expand("pair_corr/{phenocode1}.{phenocode2}.log", lambda x, y: [(x[i], y[j]) for i, j in itertools.combinations(range(0, len(x)), 2)], phenocode1 = get_all_phenocodes("pheno-list.json"), phenocode2 = get_all_phenocodes("pheno-list.json"))
   output: "result/ALL.RG.txt"
   run:
      merged = []
      for f in input:
         corr = read_corr(f)
         if not corr: continue
         merged.append(corr)
      write_merged(merged, output[0])

 
rule corr:
   input: "munged/{phenocode1}.sumstats.gz", "munged/{phenocode2}.sumstats.gz"
   output: "pair_corr/{phenocode1}.{phenocode2}.log"
   run:
      ldsc_exec = os.path.join(config['LDSC'], 'ldsc.py')
      ldsc_scores = config['LDSC_scores']
      output_prefix = os.path.join("pair_corr", "{}.{}".format(wildcards.phenocode1, wildcards.phenocode2))
      shell("{ldsc_exec} --rg {input[0]},{input[1]} --ref-ld-chr {ldsc_scores} --w-ld-chr {ldsc_scores} --out {output_prefix}")


rule munge:
   input: "traits/{phenocode}.json"
   output: "munged/{phenocode}.sumstats.gz"
   run:
      with open(input[0]) as f_in:
         trait = json.load(f_in)
         trait_phenocode = trait['phenocode']
         trait_assoc_file = trait['assoc_files'][0]
         trait_n_cases = trait['num_cases']
         trait_n_controls = trait['num_controls']
         munge_exec = os.path.join(config['LDSC'], 'munge_sumstats.py')
         snp_list = config['LDSC_snplist']
         pval_column = config['columns']['pvalue']
         effect_column = config['columns']['effect']
         snp_column = config['columns']['snp']
         effect_allele_column = config['columns']['effect_allele']
         other_allele_column = config['columns']['other_allele']
         no_effect_value = config['no_effect']
         output_prefix = os.path.join("munged", wildcards.phenocode)
         shell("{munge_exec} --sumstats {trait_assoc_file} --merge-alleles {snp_list} --N-cas {trait_n_cases} --N-con {trait_n_controls} --p {pval_column} --signed-sumstats {effect_column},{no_effect_value} --snp {snp_column} --a1 {effect_allele_column} --a2 {other_allele_column} --out {output_prefix}")


rule split:
   input: "pheno-list.json"
   output: dynamic("traits/{phenocode}.json")
   run:
      for f in input:
         with open(f) as f_in:
            data = json.load(f_in)
            for trait in data:
               with open("traits/{}.json".format(trait['phenocode']), "w") as f_out:
                  json.dump(trait, f_out)

