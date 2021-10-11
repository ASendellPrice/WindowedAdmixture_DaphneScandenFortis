#!/usr/bin/env Rscript

#Load required packages
library(pophelper)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggthemes)

#Create a list of "Q" files
Q.files <- list.files(pattern = ".Q")

#Use alingK to fix possible label switching
xlist <- alignK(readQ(Q.files), type = "within")

#Load window list
windows <- read.delim("../../window.ranges.txt", header = FALSE, sep = ":")
colnames(windows) <- c("start","end","mid")

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

        x <- paste0("window.",windows$start[WINDOW],"-",windows$end[WINDOW],".2.Q")
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


################################################################################
# Plot generation
################################################################################

MeanClusterAssigment <- read.delim("../../summarise_plot/species_cluster_assignment.txt", header = FALSE)
colnames(MeanClusterAssigment) <- c("species","cluster1_mean","cluster2_mean","mid")

MeanClusterAssigment$cluster <- NA
for (ROW in 1:nrow(MeanClusterAssigment)){
  if (MeanClusterAssigment$cluster1_mean[ROW] >= 0.5){
    MeanClusterAssigment$cluster[ROW] <- 1
  } else {
    MeanClusterAssigment$cluster[ROW] <- 2
  }
}

MeanClusterAssigment$mid <- MeanClusterAssigment$mid / 1000000


P <- ggplot(MeanClusterAssigment, aes(mid, species, fill= cluster1_mean)) +
  geom_tile()+
  #scale_fill_gradient2(low = "red",mid = "white",high = "blue",midpoint = 0.5, limits=c(0, 1)) +
  scale_fill_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red",limits=c(0, 1)) +
  #scale_fill_distiller(palette = "Spectral", limits=c(0, 1)) +
  xlab(label = "Chr1A (Mb)") +
  labs(fill="% scandens ancestry")



MeanCluster1Assigment <- MeanClusterAssigment %>%
  select(species, cluster1_mean, mid)

MeanCluster2Assigment <- MeanClusterAssigment %>%
  select(species, cluster2_mean, mid)

colnames(MeanCluster1Assigment) <- c("species","proportion","mid")
colnames(MeanCluster2Assigment) <- c("species","proportion","mid")
MeanCluster1Assigment$cluster <- "cluster1"
MeanCluster2Assigment$cluster <- "cluster2"

comb <- rbind(MeanCluster1Assigment, MeanCluster2Assigment)

P2 <- ggplot(comb, aes(fill=cluster, y=proportion, x=mid)) +
  geom_bar(position="stack", stat="identity") +
  geom_hline(yintercept = 0.5, colour="grey25", linetype="dashed") +
  scale_x_continuous(breaks = seq(0, max(comb$mid), by = 5))


pdf("../../summarise_plot/plot.pdf", width=12, height = 8)
P2 + facet_grid(species~.) + theme_few() + xlab(label = "Chr1A (Mb)") + ylab(label = "Admixture proportion")
dev.off()
