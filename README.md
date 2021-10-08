# WindowedAdmixture_DaphneScandenFortis

The following pipeline can be used to preform Admixture analyses in sliding windows. These scripts have been developed for Darwin's finches

Simple mode:

Projection mode:

Projection 2 timeframe:


Conda environment:
```
conda create -n WindowedAdmix
conda activate WindowedAdmix
conda install -c bioconda admixture bcftools vcftools plink
```


```
source scripts/WindowedAdmix.sh \
chr1A Input.vcf \
1973-1983_scandens_fortis.txt 1984-2012_scandens_fortis.txt \
500000 500000
```
