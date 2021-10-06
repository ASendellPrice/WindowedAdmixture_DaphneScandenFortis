#!/bin/bash

#Load conda environment
conda activate WindowedAdmix

#Get variable names from sbatch submission
Chrom=${1}
VCF=${2}
SampleInfo1=${3}
SampleInfo2=${4}
WindowSize=${5}
StepSize=${6}

#Save paramaters to log file
#echo `date '+%d/%m/%Y_%H:%M:%S'` > param.log
#echo "Chrom = " $Chrom >> param.log

#Create directory for focal chrom and move into it
mkdir ${CHROM}
cd ${CHROM}

#Create sample lists from info file
cat ../${SampleInfo1} | cut -f 1 > group1.list
cat ../${SampleInfo2} | cut -f 1 > group2.list

#Create focal chromosome VCF for each group and compress
vcftools --gzvcf $VCF --chr $CHROM --keep group1.list --out ${CHROM}.group1 --recode
vcftools --gzvcf $VCF --chr $CHROM --keep group2.list --out ${CHROM}.group2 --recode
bgzip ${CHROM}.group1.recode.vcf
bgzip ${CHROM}.group2.recode.vcf

#Extract sample order from VCF file
bcftools query -l ${CHROM}.group1.recode.vcf.gz > group1.order.txt
bcftools query -l ${CHROM}.group2.recode.vcf.gz > group2.order.txt

#Find last postion in VCFs
bcftools query -f '%POS\n' ${CHROM}.group1.recode.vcf.gz | tail -n 1 > length.txt

#Run R script "getWindowRanges.R" to calculate window start and end positions
Rscript ../scripts/getWindowRanges.R $WindowSize $StepSize

#Make directories for files we will be outputting
mkdir admixture_input
mkdir admixture_input/group1
mkdir admixture_input/group2
mkdir admixture_output
mkdir admixture_output/group1
mkdir admixture_output/group1/logs
mkdir admixture_output/group2
mkdir admixture_output/group2/logs

#For each window do the following . . .
for LINE in $(cat window.ranges.txt)
do

    #Extract start, end and mid points
    START=$(echo $LINE | cut -d ":" -f 1)
    END=$(echo $LINE | cut -d ":" -f 2)
    MID=$(echo $LINE | cut -d ":" -f 3)

    #Subset VCFs keeping only SNPs within window range (START to END)
    vcftools --gzvcf ${CHROM}.group1.recode.vcf.gz --chr $CHROM --from-bp $START --to-bp $END --out admixture_input/group1/window.${START}-${END} --recode
    vcftools --gzvcf ${CHROM}.group2.recode.vcf.gz --chr $CHROM --from-bp $START --to-bp $END --out admixture_input/group2/window.${START}-${END} --recode

    #Convert VCFs to plink bed, bim, fam required by admixture
    plink --vcf admixture_input/group1/window.${START}-${END}.recode.vcf --const-fid --autosome-num 30 --make-bed --out admixture_input/group1/window.${START}-${END} --allow-extra-chr
    awk '{$1=0;print $0}'  admixture_input/group1/window.${START}-${END}.bim >  admixture_input/group1/window.${START}-${END}.bim.tmp
    mv admixture_input/group1/window.${START}-${END}.bim.tmp admixture_input/group1/window.${START}-${END}.bim

    plink --vcf admixture_input/group2/window.${START}-${END}.recode.vcf --const-fid --autosome-num 30 --make-bed --out admixture_input/group2/window.${START}-${END} --allow-extra-chr
    awk '{$1=0;print $0}'  admixture_input/group2/window.${START}-${END}.bim >  admixture_input/group2/window.${START}-${END}.bim.tmp
    mv admixture_input/group2/window.${START}-${END}.bim.tmp admixture_input/group2/window.${START}-${END}.bim

    #Run admixture
    #First run with group 1 (to train admixture)
    cd admixture_output/group1
    admixture --cv ../../admixture_input/group1/window.${START}-${END}.bed 2 > logs/window.${START}-${END}.log
    #Now run admixture again for group 2 (using allele frequencies learnt from group 1)
    cd ../group2/
    cp ../group1/window.${START}-${END}.2.P window.${START}-${END}.2.P.in
    admixture -P ../../admixture_input/group2/window.${START}-${END}.bed 2 > logs/window.${START}-${END}.log
    cd ../../

    #Remove windowed VCF files and logs as these are not needed
    rm admixture_input/*/window.${START}-${END}.recode.vcf
    rm admixture_input/*/window.${START}-${END}.log
done

#Make directory where admixture results will be summarised and plotted
#mkdir ../summarise_plot
#cd admixture_output/AllSamples

#Run CLUMPP to fix between window label switching
#Prepare CLUMPP input folder
#Rscript ../../../scripts/CLUMPP_prep.R
