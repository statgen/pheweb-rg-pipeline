# pheweb-rg-pipeline
Genetic correlation calculation pipeline via summary statistics for PheWeb

Genetic correlation on the observed scale (i.e. not liability scale)

## Required tools
- Snakemake
- LDSC (https://github.com/bulik/ldsc)

## Required data:
- Summary statistics are derived from association analyses run in primarily European-ancestry samples. 
- Summary statistics contain separate columns for effect size, p-value, effect-allele (a1), the non-effect allele (a2)
- N cases and N controls (for binary traits) are provided
- N >3K 
- From LD score regression (LDSC): 
    - w_hm3.snplist #HapMap SNPs to extract for LDSC analyses 
    - eur_w_ld_chr/ #pre-computed LD scores using 1000G Eur

## How to run

### Input

The input file format `pheno-list.json` is the same as in PheWEB data import pipeline.
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

Further details on how to create input file you can find at https://github.com/statgen/pheweb.

### Configuration

Before running the pipeline you may need to change your `config.json` file:
- Specify path to the directory with LDSC tool in the `LDSC` field.
- Specify path to HapMap SNP list (i.e. w_hm3.snplist) in `LDSC_snplist` field.
- Specify path to LD scores directory (i.e. eur_w_ld_chr/) in `LDSC_scores` field.
- Provide column names inside `columns` field collection.
- Set `no_effect` to 0 if analyzing regression coefficients and 1 if analyzing odds-rations.

### Run locally

snakemake -j [number of cores]

### Output

Your final output is in `result/ALL.RG.txt` file.

Pipeline creates directories:
- `traits`: single input JSON file for each trait
- `munged`: re-formatted input files and logs by LDSC 
- `pair_corr`: LDSC genetic correlation files for each unique traits pair
- `result`: final merged result
