#!/usr/bin/env nextflow

import groovy.json.JsonSlurper

phenotypes_json = (new JsonSlurper()).parseText(file('pheno-list.json').text)

munge_params = []
phenotypes_json.eachWithIndex{ val, num -> munge_params.add( [num, val.phenocode.replaceAll('\\.', '__'), val.assoc_files[0], val.num_cases, val.num_controls] ) } 


munge_params = Channel.from(munge_params)

munge_exec = params.LDSC + "/munge_sumstats.py"
ldsc_exec = params.LDSC + "/ldsc.py"

process munge {

   label "small_mem"

   input:
   set val(num), val(phenocode), val(filename), val(n_cases), val(n_controls) from munge_params

   output:
   set val(num), val(phenocode), file("${phenocode}.sumstats.gz") into munged

   """
   ${munge_exec} --chunksize 100000 --sumstats ${filename} --merge-alleles ${params.LDSC_snplist} --N-cas ${n_cases} --N-con ${n_controls} --p ${params.columns.pvalue} --signed-sumstats ${params.columns.effect},${params.no_effect} --snp ${params.columns.snp} --a1 ${params.columns.effect_allele} --a2 ${params.columns.other_allele} --out ${phenocode}
   """
 
}


munged_a = Channel.create()
munged_b = Channel.create()
munged.into(munged_a, munged_b)


process pair_corr {

   label = "small_mem"

   input:
   set val(num1), val(phenocode1), file(munged1), val(num2), val(phenocode2), file(munged2) from munged_a.combine(munged_b).filter{ it[0] < it[3] }

   output:
   file "${phenocode1}.${phenocode2}.log" into pair_corr

   """
   ${ldsc_exec} --rg ${munged1},${munged2} --ref-ld-chr ${params.LDSC_scores} --w-ld-chr ${params.LDSC_scores} --out ${phenocode1}.${phenocode2}
   """

}


process merge {

   label "big_mem"

   publishDir 'results'

   input:
   val files from pair_corr.collect{ "'" + it + "'" }

   output:
   file "ALL.RG.txt" into merged

   """
   #!/usr/bin/env python
   import re
   import os
   def read_pair_corr(filename):
      with open(filename, 'r') as f_in:
         while True:
            line = f_in.readline()
            if not line: break
            if line.startswith('Summary of Genetic Correlation Results'): break
         header = f_in.readline()
         if not header: return {}
         header = header.rstrip().split()
         record = f_in.readline()
         if not record: return {}
         record = dict(zip(header, record.rstrip().split()))
         for p in ['p1', 'p2']:
            record[p] = re.sub(r'\\.sumstats\\.gz\$', '', os.path.split(record[p])[-1]).replace('__', '.')
         return { 'header': header, 'record': record }
   def write_merged(data, filename):
      with open(filename, 'w') as f_out:
         header = data[0]['header']
         f_out.write('{}\\n'.format('\\t'.join(header)))
         for d in data:
            f_out.write('{}\\n'.format('\\t'.join([d['record'][h] for h in header])))
   merged = []
   for f in $files:
      corr = read_pair_corr(f)
      if not corr: continue
      merged.append(corr)
   write_merged(merged, "ALL.RG.txt")
   """
}

