#!/bin/bash

#Load conda environment
conda activate WindowedAdmix

#Get variable names from sbatch submission
CHROM=${1}
VCF=${2}
SampleInfo1=${3}
SampleInfo2=${4}
WindowSize=${5}
StepSize=${6}
OutName=${7}

#Create directory for focal chrom and move into it
mkdir ${CHROM}_${OutName}
cd ${CHROM}_${OutName}

#Create sample lists from info file
cat ../${SampleInfo1} | cut -f 1 > trainAdmix.list
cat ../${SampleInfo2} | cut -f 1 > projAdmix.list
cp ../${SampleInfo2} > sample_info.tmp

#Create focal chromosome VCF for each group and compress
vcftools --gzvcf ../$VCF --chr $CHROM --keep trainAdmix.list --out ${CHROM}.trainAdmix --recode
vcftools --gzvcf ../$VCF --chr $CHROM --keep projAdmix.list --out ${CHROM}.projAdmix --recode
bgzip ${CHROM}.trainAdmix.recode.vcf
bgzip ${CHROM}.projAdmix.recode.vcf

#Extract sample order from VCF file
bcftools query -l ${CHROM}.projAdmix.recode.vcf.gz > sample_order.tmp

#Find last postion in VCFs
bcftools query -f '%POS\n' ${CHROM}.trainAdmix.recode.vcf.gz | tail -n 1 > length.txt

#Run R script "getWindowRanges.R" to calculate window start and end positions
Rscript ../scripts/getWindowRanges.R $WindowSize $StepSize

#Make directories for files we will be outputting
mkdir admixture_input
mkdir admixture_input/trainAdmix
mkdir admixture_input/projAdmix
mkdir admixture_output
mkdir admixture_output/trainAdmix
mkdir admixture_output/trainAdmix/logs
mkdir admixture_output/projAdmix
mkdir admixture_output/projAdmix/logs

#For each window do the following . . .
for LINE in $(cat window.ranges.txt)
do

    #Extract start, end and mid points
    START=$(echo $LINE | cut -d ":" -f 1)
    END=$(echo $LINE | cut -d ":" -f 2)
    MID=$(echo $LINE | cut -d ":" -f 3)

    #Subset VCFs keeping only SNPs within window range (START to END)
    vcftools --gzvcf ${CHROM}.trainAdmix.recode.vcf.gz --chr $CHROM --from-bp $START --to-bp $END --out admixture_input/trainAdmix/window.${START}-${END} --recode
    vcftools --gzvcf ${CHROM}.projAdmix.recode.vcf.gz --chr $CHROM --from-bp $START --to-bp $END --out admixture_input/projAdmix/window.${START}-${END} --recode

    #Convert VCFs to plink bed, bim, fam required by admixture
    plink --vcf admixture_input/trainAdmix/window.${START}-${END}.recode.vcf --const-fid --autosome-num 30 --make-bed --out admixture_input/trainAdmix/window.${START}-${END} --allow-extra-chr
    awk '{$1=0;print $0}'  admixture_input/trainAdmix/window.${START}-${END}.bim >  admixture_input/trainAdmix/window.${START}-${END}.bim.tmp
    mv admixture_input/trainAdmix/window.${START}-${END}.bim.tmp admixture_input/trainAdmix/window.${START}-${END}.bim

    plink --vcf admixture_input/projAdmix/window.${START}-${END}.recode.vcf --const-fid --autosome-num 30 --make-bed --out admixture_input/projAdmix/window.${START}-${END} --allow-extra-chr
    awk '{$1=0;print $0}'  admixture_input/projAdmix/window.${START}-${END}.bim >  admixture_input/projAdmix/window.${START}-${END}.bim.tmp
    mv admixture_input/projAdmix/window.${START}-${END}.bim.tmp admixture_input/projAdmix/window.${START}-${END}.bim

    #Run admixture
    #First run with group 1 (to train admixture)
    cd admixture_output/trainAdmix
    admixture --cv ../../admixture_input/trainAdmix/window.${START}-${END}.bed 2 > logs/window.${START}-${END}.log
    #Now run admixture again for group 2 (using allele frequencies learnt from group 1)
    cd ../projAdmix/
    cp ../trainAdmix/window.${START}-${END}.2.P window.${START}-${END}.2.P.in
    admixture -P ../../admixture_input/projAdmix/window.${START}-${END}.bed 2 > logs/window.${START}-${END}.log
    cd ../../

    #Remove windowed VCF files and logs as these are not needed
    rm admixture_input/*/window.${START}-${END}.recode.vcf
    rm admixture_input/*/window.${START}-${END}.log
done

#For each species (fortis/scandens) calculate the average admixture proportions across samples for each window
mkdir ../summarise_plot
cd admixture_output/projAdmix
Rscript ../../../scripts/summarise_plot.R
