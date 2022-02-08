#!/usr/bin/Rscript
# generate a PNG plot using the pipeline read retention numbers

library(methods)
library(optparse)
library(here)

inputFile <- './Hs-UMI5RACE-IgG.accounting.csv'
outputFile <- './dataRetentionPlot1.png'
debug <- FALSE


# Set up OptionParser
option_list <- list(
    make_option(c("--inputFile"), default=inputFile, help="CSV input filename [default \"%default\"]"),
    make_option(c("--outputFile"), default=outputFile, help="PNG output filename [default \"%default\"]"),
    make_option(c("--debug"), action="store_true", default=debug, help="Print out additional debugging information [default %default]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list=option_list))


df <- read.csv(file=opt$inputFile,skip = 1,header = F)
colnames(df) <- c("stage", "count")
df$stage <- factor(df$stage, levels = df$stage)
df$proportion <- df$count/df$count[1]

png(filename = here(opt$outputFile), type='cairo',
    width = 800, height = 800, pointsize = 24)

plot(NULL, xlim=c(1,nrow(df)), ylim=c(0,1), xaxt = "n", yaxt = "n", xlab='', ylab='')
title(main=paste0("Representing ",df$count[1], " reads"), xlab='stage', ylab='proportion')

lines(df$stage, df$proportion, type = 'S', lwd=3)

### Ryan Moore https://www.tenderisthebyte.com/blog/2019/04/25/rotating-axis-labels-in-r/
axis(side = 1, labels = FALSE, tck=0)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0))

text(x = 1:length(df$stage),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3]-0.03,
     ## Use names from the data list.
     labels = as.vector(df$stage),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 45 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965
)
dev.off()
