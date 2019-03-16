#!/usr/bin/env nextflow

import groovy.json.JsonSlurper

phenotypes_json = (new JsonSlurper()).parseText(file('pheno-list.json').text)

format_params = []
phenotypes_json.eachWithIndex{ val, num -> format_params.add( [num, val.phenocode.replaceAll('\\.', '__'), val.assoc_files[0]] ) } 

bim = Channel.fromPath(params.ref_panel + "/*.bim").collect()
bed = Channel.fromPath(params.ref_panel + "/*.bed").collect()
fam = Channel.fromPath(params.ref_panel + "/*.fam").collect()

ldak_exec = params.LDAK


process mhc {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 3

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
   maxRetries 3

   input:
   set val(num), val(phenocode), val(filename) from format_params
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


process intersection {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 3

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
   maxForks 1
  
   input:
   set file(exclude1), file(exclude2) from formatted2unique_a.combine(formatted2unique_b).filter{ it[0] < it[2] }.map { [it[1], it[3]] } 

   output:
   set file(exclude1), file(exclude2), stdout into unique 

   """
   printf `cat ${exclude1} ${exclude2} | sort | uniq | md5sum | awk '{print \$1;}'`
   """

}


process tagging {
   
   label "high_cpu"
   errorStrategy "retry"
   maxRetries 3

   input:
   file intersected from intersected
   set file(exclude1), file(exclude2), val(mdsum) from unique.unique { it[2] }
   file bed from bed
   file bim from bim
   file fam from fam

   output:
   file "sumldak_${mdsum}.tagging" into tagged

   """
   cat ${exclude1} ${exclude2} | sort | uniq > combined.exclude
   for bim_file in `find . -name "*.bim"`; do
      for chr in `cut -f1 \${bim_file} | uniq`; do
         (${ldak_exec} --cut-weights weights_${mdsum}_\${chr} --bfile \${bim_file%*.bim} --extract ${intersected} --exclude combined.exclude --chr \${chr} &&\
         ${ldak_exec} --calc-weights-all weights_${mdsum}_\${chr} --bfile \${bim_file%*.bim} --extract ${intersected} --exclude combined.exclude --chr \${chr} &&\
         ${ldak_exec} --calc-tagging sumldak_${mdsum}_\${chr} --bfile \${bim_file%*.bim} --weights weights_${mdsum}_\${chr}/weights.short --power -0.25 --extract ${intersected} --exclude combined.exclude --window-kb 1000 --chr \${chr} ) &
      done
   done 
   wait
   find . -maxdepth 1 -name "sumldak_${mdsum}_*.tagging" -printf "%f\n" | sort -V > list.txt
   ${ldak_exec} --join-tagging sumldak_${mdsum} --taglist list.txt
   """

}


process pair_corr {

   label "small_mem"
   errorStrategy "retry"
   maxRetries 3

   input:
   set val(phenocode1), file(stats1), file(exclude1), val(phenocode2), file(stats2), file(exclude2) from formatted2pair_corr_a.combine(formatted2pair_corr_b).filter{ it[0] < it[4] }.map { [it[1], it[2], it[3], it[5], it[6], it[7]] }
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
   maxRetries 3

   publishDir "results"

   input:
   file files from pair_corr.collectFile() { item -> ["list.txt", "${item}\n"] } 

   output:
   file "ALL.RG.txt" into merged

   """
   merge.py -i ${files} -o ALL.RG.txt
   """
}

