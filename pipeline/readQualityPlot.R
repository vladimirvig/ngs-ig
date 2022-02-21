#!/usr/bin/Rscript

library(methods)
library(dada2)

inputFile <- ''
outputFile <- ''
debug <- FALSE

# Set up OptionParser
option_list <- list(
  optparse::make_option(c("--inputFile"), default=inputFile, help="Absolute path of the input FASTQ file [default \"%default\"]"),
  optparse::make_option(c("--outputFile"), default=outputFile, help="Absolute path of the output PDF file [default \"%default\"]"),
  optparse::make_option(c("--debug"), action="store_true", default=debug, help="Print out additional debugging information [default %default]")
)

# Parse arguments
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

pdf(opt$outputFile)
plotQualityProfile(opt$inputFile)
dev.off()
