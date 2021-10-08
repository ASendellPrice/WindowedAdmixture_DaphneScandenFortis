# WindowedAdmixture_DaphneScandenFortis

## Introduction
The following pipeline can be used to preform Admixture projection analyses in sliding windows. These were developed for the Darwin's finches project to assess population structuring within a dataset with un-even sampling between focal species. The pipeline is developed for use with two species/population groupings: here we use Daphne scandens and Daphne fortis.

The main script "WindowedAdmix.sh" consists of two main components:

1. An initial admixture run using an equal number of samples per species.
2. A second admixture run on the full target samples using the allele frequencies learned during step 1.

## Environement setup
To run the pipeline, the following conda environment needs to be set up:
```
conda create -n WindowedAdmix
conda activate WindowedAdmix
conda install -c bioconda admixture bcftools vcftools plink
```

Note: R will also need to be installed, along with the following packages:
- dplyr
- pophelper
- ggplot2
- ggthemes

## Running the pipeline
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
