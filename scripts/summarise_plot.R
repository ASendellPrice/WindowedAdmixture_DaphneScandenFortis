#!/usr/bin/env Rscript

#Load required packages
library(pophelper)
library(dplyr)

#Create a list of "Q" files
Q.files <- list.files(pattern = ".Q")

#Use alingK to fix possible label switching
xlist <- alignK(readQ(Q.files), type = "within")

#Load window list
windows <- read.delim("../../window.ranges.txt", header = FALSE, sep = " ")
colnames(windows) <- c("start","end","mid")

#Remove windows with zero SNPs
#windows <- filter(windows, sites > 0)

#Load sample ID and species info
sample_info <- read.delim("../../sample_info.tmp", header = FALSE)
colnames(sample_info) <- c("sampleID","location","species")

#Load sample order in VCF
sample_order <- read.delim("../../sample_order.tmp", header = FALSE)
colnames(sample_order) <- c("sampleID")

#Join sample order and sample info so that they are in the correct order
combined <- right_join(sample_order, sample_info)

#For each window do the following . . .
for (WINDOW in 1:nrow(windows)){

        x <- paste0("window.",windows$start[WINDOW],"-",windows$end[WINDOW],".AllSamples.2.Q")
        window.Q <- cbind(combined, as.data.frame(xlist[[x]][[1]]), as.data.frame(xlist[[x]][[2]]))
        colnames(window.Q) <- c("sampleID","location","species","cluster1","cluster2")
        window.Q$mid <- windows$mid[WINDOW]


        MeanClusterAssigment <- window.Q %>%
        group_by(species) %>%
        summarise_at(vars(cluster1, cluster2), list(mean = mean))
        MeanClusterAssigment$mid <- windows$mid[WINDOW]

        write.table(MeanClusterAssigment, "../../summarise_plot/species_cluster_assignment.txt",
                sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
}

MeanClusterAssigment <- read.delim("species_cluster_assignment.txt", header = FALSE)
colnames(MeanClusterAssigment) <- c("species","cluster1_mean","cluster2_mean","mid")
