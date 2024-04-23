# fairy-test

Scripts for plotting the benchmarking results for fairy. 

## Requirements

1. seaborn
2. matplotlib
3. numpy
4. scipy
5. pandas

## Snakemake pipeline for generating results

Results were generated using the `bin_mags.smk` pipeline. This snakemake file was customized depending on the assembler/dataset, but all of the rules were the same across all datasets. 

## Figure 1

To regenerate results for Figure 1, run scripts/pearson.py and scripts/timeplot.py scripts.

```sh
python scripts/pearson.py
python scripts/timeplot.py
```

### Figure 3A/B

To regenerate results for Figure 2B, run scripts/binnercomp.py.

```sh
python scripts/binnercomp.py
```

### Figure 2

To regenerate Figure 2A, run any combination of 
```sh
python scripts/fig2.py chicken_new biofilm_new ...
```

`*_new` indicates the checkm2 results for MetaBAT2. `*_maxbin` are checkm2 results for MaxBin2, and etc. 

## Figure 4

Check out the `notebooks/contamination.ipynb` jupyter notebook. You will need jupyter notebook installed. 

