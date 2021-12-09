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

# 3.load database annotation
#Load biomaRt library and you should be connected to internet
library(biomaRt) 

# 4. Choose mart to use
#Check the list mart available and choose which mart we are interested in our project

listMarts()


#We want to use mart "ENSEMBL_MART_ENSEMBL in our project.
ensembl = useMart("ENSEMBL_MART_ENSEMBL") #this line tell BioMart to connect to this specific mart

#In each specific mart under BioMart there are different species by each dataset. 
#The list of datasets with descriptions and version.
listDatasets(ensembl, verbose = TRUE)


#Load the dataset of choice and collect gene names from homo sapiens dataset
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
