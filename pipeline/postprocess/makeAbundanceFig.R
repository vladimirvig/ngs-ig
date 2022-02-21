#!/usr/bin/Rscript
# generate a PDF plot for binned population analysis for UMI-tagged NGS data

library(methods)
library(ggplot2)

f <- file ("stdin")
open(f)
data <- readLines(f)
close (f)

if(is.na(as.numeric(data[1]))) { stop("Error: input is not a column of integers; please check!") }

df <- as.data.frame(as.numeric(data))
df1<- data.frame(summary(cut(df[,1], c(0,2^(0:17)), labels = 2^(0:17))))

total <- sum(df)

colnames(df1) <- c("N")
df1$bins <- as.numeric(rownames(df1))

df1$abundance <- df1$N * df1$bins
df1$proportion <- df1$abundance/sum(df1$abundance)

df1 <- df1[df1$bins>1,]

p <- ggplot(df1, aes(x=bins))

p <- p + geom_bar(aes(y=N/max(N)/2), stat = "identity", fill="blue", alpha = 0.5)

p <- p + geom_line(aes(y=proportion, color = "proportion")) + scale_x_log10(name = "MIG size, reads", expand=c(0,0), limits=c(1, max(df1$bins)*2), breaks = c(1:10,100,1000,10000,100000), labels = c("1", rep("", 8), 10, 100, 1000, 10000, 100000), oob=scales::rescale_none)  + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle (label = paste0("Representing ", total, " reads"))

p <- p + scale_y_continuous(name = "proportion", sec.axis = sec_axis(~.*2*max(df1$N), name = "abundance")) + theme(legend.position = c(0.5, 0.9), axis.text.y.right = element_text(color = "blue"), axis.line.y.right = element_line(color = "blue"), axis.title.y.right=element_text(color = "blue")) + scale_colour_manual(values = c("black", "blue")) + labs (colour = "Parameter")

pdf("output.pdf")
p
dev.off()
