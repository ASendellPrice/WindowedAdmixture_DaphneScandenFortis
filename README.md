# WindowedAdmixture_DaphneScandenFortis
Pipeline for running window based admixture analyses for finches

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
