# WindowedAdmixture_DaphneScandenFortis

The following pipeline can be used to preform Admixture analyses in sliding windows. These scripts have been developed for Darwin's finches

Simple mode:

Projection mode:

Projection 2 timeframe:


Conda environment: ADD R
```
conda create -n WindowedAdmix
conda activate WindowedAdmix
conda install -c bioconda admixture bcftools vcftools plink
```


To perform windowed admixture analysis submit script along with the following parameters:
1. name of focal chromosome (as appears in VCF)
2. input VCF file
3. training sample list with equal numbers of species A and species B
4. target sample info specifying: sampleID  Island  Species
5. window size in bp
6. step sixe in bp
7. suffix for output directory name

```
source scripts/WindowedAdmix.sh \
chr1A Input.vcf \
31S_31F_to1983.txt 1973-1983_scandens_fortis.txt \
500000 500000 \
Daphne_Sc_Fo_to1983
```
