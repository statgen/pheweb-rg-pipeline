# pheweb-rg-pipeline

Pipeline for calculating genetic correlations via summary statistics between >1,000 phenotypes in PheWeb.
Genetic correlation is on the observed scale (i.e. not liability scale).

Pipeline allows to choose the following tools:
- LDSC (https://github.com/bulik/ldsc)
- SumHer (http://dougspeed.com/sumher/)

Pipeline can be run locally or on SLURM.

## Required tools
- Python 3 (recommended)
- Nextflow (https://www.nextflow.io) 
   * can be installed as a standalone tool (https://www.nextflow.io/docs/latest/getstarted.html#installation), or
   * using Miniconda (https://anaconda.org/bioconda/nextflow)
- LDSC (https://github.com/bulik/ldsc), when computing genetic correlations via LDSC 
- SumHer (http://dougspeed.com/sumher/), when computing genetic correlations via SumHer

## Required data:
- Summary statistics are derived from association analyses run in primarily European-ancestry samples. 
- Summary statistics contain separate columns for effect size, p-value, effect-allele (a1), the non-effect allele (a2)
- N cases and N controls (for binary traits) are provided
- N >3K 
- When using LDSC: 
    - w_hm3.snplist #HapMap SNPs to extract for LDSC analyses (see https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
    - eur_w_ld_chr/ #pre-computed LD scores using 1000G Eur (see https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
- When using SumHer:
    - 1000 Genomes based reference panel in PLINK binary format (i.e. can use the 1000 Genome EUR phase 3 files provided by LDSC at https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz) 

## How to run

### Input

The input file must be named `pheno-list.json`. It has the same format as in the PheWeb data import pipeline.
```json
[
 {
  "assoc_files": ["/home/watman/ear-length.epacts.gz"],
  "phenocode": "ear-length",
  "num_cases": 10000,
  "num_controls": 100
 },
 {
  "assoc_files": ["/home/watman/eats-kimchi.autosomal.epacts.gz"],
  "phenocode": "eats-kimchi",
  "num_cases": 14000,
  "num_controls": 100
 }
]
```

The key difference are:
-  Fields `num_cases` and `num_controls` are required
-  Only one file per trait is allowed in `assoc_files` (i.e. cannot split the summary statistics by chromosome)

If you have a tab-delimited file (no header) with the following columns: full path to summary stat file, phenocode, number of cases, number of controls, use `tab2pheno-list.py -i [file.tsv]` to create `pheno-list.json`.

Further details on how to create the input file are at https://github.com/statgen/pheweb.

### Run LDSC

#### - Configuration

Before running the pipeline you may need to change your `LDSC/nextflow.config` file:
- Specify path to the directory where LDSC is installed in the `LDSC` field.
- Specify path to the HapMap SNP list (i.e. w_hm3.snplist) in the `LDSC_snplist` field.
- Specify path to the LD scores directory (i.e. eur_w_ld_chr/) in the `LDSC_scores` field.
- Provide column names inside the `columns` configuration scope.
- Set `no_effect` to 0 if analyzing regression coefficients or 1 if analyzing odds-ratios.

#### - Locally 

Inside the `LDSC/nextflow.config` file, set the number of cpus you want to use via the `cpus` parameter:
```
...
$local {
  cpus = 4
}
...
```

Place your input file `pheno-list.json` inside the directory where you want to save results (this will also be the working directory for all intermediate files). Then, in the same directory run:
```
nextflow run /path/to/LDSC.nf
```

#### - With SLURM

Inside the `LDSC/nextflow.config` file, uncomment `executor = "slurm"` line and comment `executor = "local"` line:
```
executor = "slurm"
// executor = "local"
```
Set SLURM queue name via the `queue` parameter.

Place your input file `pheno-list.json` inside the directory where you want to save results (this will also be the working directory for all intermediate files). Then, in the same directory run:
```
nextflow run /path/to/LDSC.nf
```

### Run SumHer

#### - Configuration

Before running the pipeline you may need to change your `SumHer/nextflow.config` file:
- Specify path to the LDAK executable in the `LDAK` field.
- Specify path to the directory with the referfernce panel (in Plink's bim and bam files; may be split by chromosome) in the `ref_panel` field.
- (For incremental apporach) Specify the `*.json` file with a subset of summary statistics from `pheno-list.json` in the `compute_pheno_list` field.

  Rationale:
  
  Sometimes you may need to compute genetic correlations for summary statistics files which arrive in batches. To avoid  recomputing all genetic correlations every time, our pipeline supports "incremental" approach. Suppose you had an input file `pheno-list.json` with 10 summary statistics files and you computed genetic correlations using this file (45 unique pairs). Later you got 4 other summary statistics files and created a new input file `new-pheno-list.json`. One approach is to combine 
`pheno-list.json` and `new-pheno-list.json` together and run the pipeline again. In this way you will perform unnecessary computations of genetic correlations that you computed previously (45 out of 91 correlations). However, if you specify `new-pheno-list.json` in the `compute_pheno_list` field, then the pipeline will know that you want to compute only 46 new correlations with and within your 4 new summary statistics files. This can safe significant amount of computational time when you are dealing with >1000 phenotypes.

  To combine two (or more) phenotype files use `cat-pheno-list.py -i pheno-list.json new-pheno-list.json -o pheno-list.json`. Note: if you don't want to overwrite your original `pheno-list.json` file, specify different name in `-o` option and change the corresponding name in the `all_pheno_list` field of your `SumHer/nextflow.config` file.
  
  Chunking:
  
  Alternatively, you can chunk a large `pheno-list.json` file into multiple independent non-overlapping runs.
  For example, the following command chunks `pheno-list.json` into chunks with 100 phenotypes each:
  ```
  chunk-pheno-list.py -i pheno-list.json -s 100
  ```
  It will create multiple `chunk_[1-9]+` directories. Each directory will store `pheno-list.json` and `compute-pheno-list.json` files.
  In order to run a single chunk (for example `chunk_1`):
  - set `compute_pheno_list` to `compute-pheno-list.json` in your `SumHer/nextflow.config` file
  - `cd chunk_1`
  - `nextflow run /path/to/SumHer/SumHer.nf`
  

#### - Locally 

Inside the `SumHer/nextflow.config` file, set the number of cpus you want to use via the `cpus` parameter e.g.:
```
...
$local {
  cpus = 4
}
...
```

Place your input file `pheno-list.json` inside the directory where you want to save results (this will also be the working directory for all intermediate files). Then, in the same directory run:
```
nextflow run /path/to/SumHer.nf
```

#### - With SLURM

Inside the `SumHer/nextflow.config` file:
1. uncomment `executor = "slurm"` line and comment `executor = "local"` line e.g.:
  ```
  executor = "slurm"
  // executor = "local"
  ```
2. Set SLURM queue name via the `queue` parameter.
3. Set maximal number of parallel SLURM jobs via the `queueSize` e.g.:
  ```
  $slurm {
    queueSize = 1000
  }
  ```

Place your input file `pheno-list.json` inside the directory where you want to save results (this will also be the working directory for all intermediate files). Then, in the same directory run:
```
nextflow run /path/to/LDSC.nf
```


### Output

Your final output (matrix of correlations among all the traits) is in the `workdir` directory in the `result/ALL.RG.txt` file.

The pipeline creates directories:
- `work`: `Nexflow` working directory with output files from all steps.
- `result`: directory with final merged result
