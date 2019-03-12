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
- LDSC (https://github.com/bulik/ldsc), when computing LDSC genetic correlations
- LDAK (http://dougspeed.com/sumher/), when computing SumHer genetic correlations

## Required data:
- Summary statistics are derived from association analyses run in primarily European-ancestry samples. 
- Summary statistics contain separate columns for effect size, p-value, effect-allele (a1), the non-effect allele (a2)
- N cases and N controls (for binary traits) are provided
- N >3K 
- When using LDSC: 
    - w_hm3.snplist #HapMap SNPs to extract for LDSC analyses (see https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
    - eur_w_ld_chr/ #pre-computed LD scores using 1000G Eur (see https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
- When using SumHer:
    - 1000 Genomes based reference panel from ...

## How to run

### Configuration

Before running the pipeline you may need to change your `nextflow.config` file:
- Specify path to the directory where LDSC is installed in the `LDSC` field.
- Specify path to the HapMap SNP list (i.e. w_hm3.snplist) in the `LDSC_snplist` field.
- Specify path to the LD scores directory (i.e. eur_w_ld_chr/) in the `LDSC_scores` field.
- Provide column names inside the `columns` configuration scope.
- Set `no_effect` to 0 if analyzing regression coefficients or 1 if analyzing odds-ratios.

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

#### - Locally 

Inside the `SumHer/nextflow.config` file, set the number of cpus you want to use via the `cpus` parameter:
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
[ under preparation ]


### Output

Your final output (matrix of correlations among all the traits) is in the `workdir` directory in the `result/ALL.RG.txt` file.

The pipeline creates directories:
- `work`: `Nexflow` working directory with output files from all steps.
- `result`: directory with final merged result
