# pheweb-rg-pipeline
Genetic correlation calculation pipeline via summary statistics for PheWeb

Genetic correlation on the observed scale (i.e. not liability scale)

#### Pre-requisites:
- Summary statistics are derived from association analyses run in primarily European-ancestry samples. 
- Summary statistics contain separate columns for effect size, p-value, effect-allele (a1), the non-effect allele (a2)
- N cases and N controls (for binary traits) are provided
- N >3K 
- From LD score regression (LDSC): 
    - munge_sumstats.py #LDSC file formatting script
    - w_hm3.snplist #HapMap SNPs to extract for LDSC analyses 
    - eur_w_ld_chr/ #pre-computed LD scores using 1000G Eur
    - ldsc.py #LDSC genetic correlation script
