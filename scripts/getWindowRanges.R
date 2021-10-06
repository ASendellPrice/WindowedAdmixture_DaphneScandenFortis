#!/usr/bin/env Rscript

#Keep command line arguments
args <- commandArgs()

# set chromosome size
chrom_length <- as.numeric(readLines("length.txt"))

# set window size and window jump
window_size <- args[1]
window_jump <- args[2]

# use seq to find the start points of each window
window_start <- seq(from = 1, to = chrom_length, by = window_jump)
# add the size of the window to each start point
window_stop <- window_start + window_size-1

# no windows start before the end of chromosome
sum(window_start > chrom_length)
# but some window stop positions do occur past the final point
sum(window_stop > chrom_length)

# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < chrom_length)]
window_stop <- window_stop[which(window_stop < chrom_length)]

# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop,
                      mid = window_start + (window_stop-window_start)/2)

write.table(windows, "window.ranges.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
