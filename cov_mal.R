# 1. set working directory
setwd("D:/Pdf/Codeathon Project/counts/counts1")


# path to save the outputs
outpath <- "D:/Pdf/Codeathon Project/counts/counts1/"  


# 2. Import data
#matrix tells how many reads have been mapped to a gene.
counts1 <- read.delim("counts.matrix1.renamed", sep="\t", header = TRUE, row.names = 1) #counts matrix containing the read counts.

head(counts1)

dim(counts1)

colnames(counts1)

#delete col3 (control3) 
counts1 <- subset(counts1, select = -control3)

dim(counts1)

head(counts1)


colnames(counts1)

summary(counts1)