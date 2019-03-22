#!/usr/bin/env nextflow

import groovy.json.JsonSlurper

phenotypes_json = (new JsonSlurper()).parseText(file(params.all_pheno_list).text)

compute_phenotypes_json = []
if (params.compute_pheno_list != "") {
   compute_phenotypes_json = (new JsonSlurper()).parseText(file(params.compute_pheno_list).text)
}

if (!phenotypes_json.containsAll(compute_phenotypes_json)) {
   println "The ${params.all_pheno_list} file must contain all entries from the ${params.compute_pheno_list} file (i.e. must be a superset)."
   exit 1
}

all_phenotypes = []
phenotypes_json.eachWithIndex{ val, num -> all_phenotypes.add( [num, val.phenocode.replaceAll('\\.', '__'), val.assoc_files[0]] ) } 

compute_phenotypes = []
compute_phenotypes_json.each {
   phenocode = it.phenocode.replaceAll('\\.', '__');
   num = all_phenotypes.find { it[1] == phenocode }[0]; // there will be always a match because of the previous containsAll check
   compute_phenotypes.add(num);
}

ldak_exec = params.LDAK


bim = Channel.fromPath(params.ref_panel + "/*.bim").collect()
bed = Channel.fromPath(params.ref_panel + "/*.bed").collect()
fam = Channel.fromPath(params.ref_panel + "/*.fam").collect()


process mhc {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   input:
   file bim from bim

   output:
   file "mhc.variants" into mhc

   """
   mhc.py -i ${bim} -b GRCh37 -o mhc.variants 
   """

}


process format {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   input:
   set val(num), val(phenocode), val(filename) from all_phenotypes
   file bed from bed
   file bim from bim
   file fam from fam
   file mhc from mhc

   output:
   set val(num), val(phenocode), file("${phenocode}.stats"), file("${phenocode}.nonamb"), file("${phenocode}.exclude") into formatted

   """
   format.py -i ${filename} -o ${phenocode}
   if [ -s ${phenocode}.big ]; then
      i=0
      for f in ${bed}; do
         i=\$((i+1))
         ${ldak_exec} --remove-tags ${phenocode}.\${i} --bfile \${f%*.bed} --top-preds ${phenocode}.big --window-kb 1000 --min-cor 0.1   
      done
      cat ${phenocode}.*.out > ${phenocode}.out
   else
      touch ${phenocode}.out
   fi
   cat ${mhc} ${phenocode}.out > ${phenocode}.exclude
   """
 
}


formatted2intersection = Channel.create()
formatted2unique_a = Channel.create()
formatted2unique_b = Channel.create()
formatted2pair_corr_a = Channel.create()
formatted2pair_corr_b = Channel.create()


formatted.separate(formatted2intersection, formatted2unique_a, formatted2unique_b, formatted2pair_corr_a, formatted2pair_corr_b) { 
   it -> [ it[3], [it[0], it[4]], [it[0], it[4]], [it[0], it[1], it[2], it[4]], [it[0], it[1], it[2], it[4]]  ]
}

if (compute_phenotypes) {
   unique_input = formatted2unique_a.filter { compute_phenotypes.contains(it[0]) }
      .combine(formatted2unique_b).filter{  !compute_phenotypes.contains(it[2]) ? true : it[0] < it[2] }.map { [it[1], it[3]] }
   pair_corr_input = formatted2pair_corr_a.filter { compute_phenotypes.contains(it[0]) }
      .combine(formatted2pair_corr_b).filter{ !compute_phenotypes.contains(it[4]) ? true : it[0] < it[4] }.map { [it[1], it[2], it[3], it[5], it[6], it[7]] }
} else {
   unique_input = formatted2unique_a.combine(formatted2unique_b).filter{ it[0] < it[2] }.map { [it[1], it[3]] }
   pair_corr_input = formatted2pair_corr_a.combine(formatted2pair_corr_b).filter{ it[0] < it[4] }.map { [it[1], it[2], it[3], it[5], it[6], it[7]] }
}


process intersection {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   input:
   file nonambs from formatted2intersection.collect()

   output:
   file ("intersection.nonamb") into intersected

   """
   intersect.py -i ${nonambs} -o "intersection.nonamb"
   """

}


process unique {

   executor 'local'
  
   input:
   set file(exclude1), file(exclude2) from unique_input

   output:
   set file(exclude1), file(exclude2), stdout into unique 

   """
   printf `cat ${exclude1} ${exclude2} | sort | uniq | md5sum | awk '{print \$1;}'`
   """

}


process tagging_chr {
   
   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   input:
   file intersected from intersected
   set file(exclude1), file(exclude2), val(mdsum) from unique.unique { it[2] }
   file bed from bed
   file bim from bim
   file fam from fam
   each chr from Channel.from(1..22)

   output:
   set val(mdsum), file("sumldak_${mdsum}_${chr}.tagging") into tagged_chr

   """
   cat ${exclude1} ${exclude2} | sort | uniq > combined.exclude
   bim_file=`grep -l "^${chr}\\s" *.bim`
   [ -z \$bim_file ] && echo "BIM file for chromosome ${chr} was not found!" && exit 1
   ${ldak_exec} --cut-weights weights_${mdsum}_${chr} --bfile \${bim_file%*.bim} --extract ${intersected} --exclude combined.exclude --chr ${chr}
   ${ldak_exec} --calc-weights-all weights_${mdsum}_${chr} --bfile \${bim_file%*.bim} --extract ${intersected} --exclude combined.exclude --chr ${chr}
   ${ldak_exec} --calc-tagging sumldak_${mdsum}_${chr} --bfile \${bim_file%*.bim} --weights weights_${mdsum}_${chr}/weights.short --power -0.25 --extract ${intersected} --exclude combined.exclude --window-kb 1000 --chr ${chr} 
   """
 
}


process tagging_merge {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   input:
   set val(mdsum), file(tagged_chr) from tagged_chr.groupTuple()

   output:
   file "sumldak_${mdsum}.tagging" into tagged

   """
   find . -maxdepth 1 -name "sumldak_${mdsum}_*.tagging" -printf "%f\n" | sort -V > list.txt
   ${ldak_exec} --join-tagging sumldak_${mdsum} --taglist list.txt
   """

}


process pair_corr {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   input:
   set val(phenocode1), file(stats1), file(exclude1), val(phenocode2), file(stats2), file(exclude2) from pair_corr_input
   file tagged from tagged.collect()

   output:
   file "${phenocode1}.${phenocode2}.cors" into pair_corr

   """
   mdsum=`cat ${exclude1} ${exclude2} | sort | uniq | md5sum | awk '{print \$1;}'`
   tagfile="sumldak_\${mdsum}.tagging"
   ${ldak_exec} --sum-cors ${phenocode1}.${phenocode2} --tagfile \${tagfile} --summary ${stats1} --summary2 ${stats2} --genomic-control YES --check-sums NO
   """

}


process merge {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 5

   publishDir "results"

   input:
   file files from pair_corr.collectFile() { item -> ["list.txt", "${item}\n"] } 

   output:
   file "ALL.RG.txt" into merged

   """
   merge.py -i ${files} -o ALL.RG.txt
   """
}
