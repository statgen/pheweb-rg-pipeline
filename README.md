# pheweb-rg-pipeline

Pipepline for calculating genetic correlations via summary statistics between >1,000 phenotypes in PheWeb.
Genetic correlation is on the observed scale (i.e. not liability scale)

## Required tools
- Snakemake
- LDSC (https://github.com/bulik/ldsc)

## Required data:
- Summary statistics are derived from association analyses run in primarily European-ancestry samples. 
- Summary statistics contain separate columns for effect size, p-value, effect-allele (a1), the non-effect allele (a2)
- N cases and N controls (for binary traits) are provided
- N >3K 
- From LD score regression (LDSC): 
    - w_hm3.snplist #HapMap SNPs to extract for LDSC analyses (see https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
    - eur_w_ld_chr/ #pre-computed LD scores using 1000G Eur (see https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)

## How to run

### Configuration

Before running the pipeline you may need to change your `config.json` file:
- Set your working directory in the `workdir` field. Pipeline will write all intermediate and final results to this working directory. Default working directory is `.` (i.e. current directory with the `Snakemake` file).
- Specify path to the directory with LDSC tool in the `LDSC` field.
- Specify path to HapMap SNP list (i.e. w_hm3.snplist) in `LDSC_snplist` field.
- Specify path to LD scores directory (i.e. eur_w_ld_chr/) in `LDSC_scores` field.
- Provide column names inside `columns` field collection.
- Set `no_effect` to 0 if analyzing regression coefficients and 1 if analyzing odds-rations.

### Input

The input file must be named `pheno-list.json` and must be placed to your `workdir` (see Configuration section). It has the same format as in PheWEB data import pipeline.
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
-  Only one file is allowed in `assoc_files`

If you have a tab-delimited file (no header) with the following columns: full path to summary stat file, phenocode, number of cases, number of controls, use `tab2pheno-list.py -i [file.tsv]` to create `pheno-list.json`.

Further details on how to create input file you can find at https://github.com/statgen/pheweb.

### Run

*Important:* if you have >1000 phenotypes, it might take a while until `snakemake` evaluates all dependenies and starts submitting computational jobs. Try to use the most recent version of `snakemake`.

#### - Locally 

```
snakemake -j [number of cores to use]
```

#### - With SLURM

(Optional) Edit `partition` option in the `cluster.SLURM.json` configuration file.

Run the following command specifying max. number of SLURM jobs (submitted at once) with the `-j` option:
```
snakemake -j [max number of jobs to submit] --cluster-config cluster.SLURM.json --cluster "sbatch -J {cluster.job-name} --mem {cluster.mem} -p {cluster.partition} -t {cluster.time} -e {cluster.error} -o {cluster.output}"
```

### Output

Your final output is in the `workdir` directory in the `result/ALL.RG.txt` file.

Pipeline creates directories:
- `traits`: single input JSON file for each trait
- `munged`: re-formatted input files and logs by LDSC 
- `pair_corr`: LDSC genetic correlation files for each unique traits pair
- `result`: final merged result
