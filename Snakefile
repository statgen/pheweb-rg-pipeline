import json
import os
import itertools


configfile: "config.json"
workdir: config.get('workdir', '.')


def get_phenocodes(json_file):
   with open(json_file, 'r') as f_in:
      traits = json.load(f_in)
      return [str(t['phenocode']).replace('.', '__') for t in traits]


def get_phenocode_buddy(phenocodes, phenocode):
   for i in range(phenocodes.index(phenocode) + 1, len(phenocodes)):
      yield phenocodes[i]


def get_munge_params(filename, phenocode):
   with open(filename) as f_in:
      traits = json.load(f_in)
      for trait in traits:
         munged_phenocode = str(trait['phenocode']).replace('.', '__')
         if phenocode == munged_phenocode:
            return [trait['assoc_files'][0], trait['num_cases'], trait['num_controls']]


phenocodes = get_phenocodes("pheno-list.json")
ldsc_exec = os.path.join(config['LDSC'], 'ldsc.py')
munge_exec = os.path.join(config['LDSC'], 'munge_sumstats.py')


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
         record[p] = re.sub(r'\.sumstats\.gz$', '', os.path.split(record[p])[-1]).replace('__', '.')
      return {'header': header, 'record': record}


def write_merged(data, filename):
   with open(filename, 'w') as f_out:
      header = data[0]['header']
      f_out.write('{}\n'.format('\t'.join(header)))
      for d in data:
         f_out.write('{}\n'.format('\t'.join([d['record'][h] for h in header])))


rule all:
   input: "result/ALL.RG.txt"


rule merge_all:
   input: expand("pair_corr/{phenocode}.RG.txt", phenocode = phenocodes)
   output: "result/ALL.RG.txt"
   run:
      with open(output[0], 'w') as f_out:
         for i, f in enumerate(input):
            with open(f, 'r') as f_in:
               header = f_in.readline()
               if not header: continue
               if i == 0: f_out.write(header)
               for line in f_in: f_out.write(line)


rule merge_single:
   input: lambda wildcards: [ "pair_corr/{0}/{0}.{1}.log".format(wildcards.phenocode, p)  for p in get_phenocode_buddy(phenocodes, wildcards.phenocode) ]
   output: "pair_corr/{phenocode}.RG.txt"
   run:
      if not input: # this phenocode has no unique buddies
         shell("touch {output}")
      else:
         merged = []
         for f in input:
            corr = read_corr(f)
            if not corr: continue
            merged.append(corr)
         write_merged(merged, output[0])

 
rule corr:
   input: "munged/{phenocode1}.sumstats.gz", "munged/{phenocode2}.sumstats.gz"
   params:
      prefix = "pair_corr/{phenocode1}/{phenocode1}.{phenocode2}"
   output: "pair_corr/{phenocode1}/{phenocode1}.{phenocode2}.log"
   conda: "ldsc_env.yml"
   shell:
      "{ldsc_exec} --rg {input[0]},{input[1]} --ref-ld-chr {config[LDSC_scores]} --w-ld-chr {config[LDSC_scores]} --out {params.prefix}"


rule munge:
   input: "pheno-list.json"
   params: 
      prefix = "munged/{phenocode}",
      munge_params = lambda wildcards, input: get_munge_params(input[0], wildcards.phenocode)
   output: "munged/{phenocode}.sumstats.gz"
   conda: "ldsc_env.yml"
   shell:
      "{munge_exec} --chunksize 100000 --sumstats {params.munge_params[0]} --merge-alleles {config[LDSC_snplist]} --N-cas {params.munge_params[1]} --N-con {params.munge_params[2]} --p {config[columns][pvalue]} --signed-sumstats {config[columns][effect]},{config[no_effect]} --snp {config[columns][snp]} --a1 {config[columns][effect_allele]} --a2 {config[columns][other_allele]} --out {params.prefix}"

