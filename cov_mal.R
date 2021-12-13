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

# 8. Extract the targeted transcripts

#Extract data from the hsapiens_gene_ensembl dataset via BioMart and get a list of transcripts.
#We have the access to the genomic data for hsapiens_gene_ensembl provided by BioMart.
#add genes names into the table
t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", 
                 "transcript_version", 
                 "ensembl_gene_id",
                 "external_gene_name",
                 "description",
                 "transcript_biotype"), mart = mart)

#Renaming table for sleuth analysis
#Combine two columns (ensembl_transcript_id and transcript_version) both separate by "." into one column called "target_id
t2g <- dplyr::mutate(t2g, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))



#Rename the table into three column such as target_id, Geneid and Gene_Name.
t2g <- dplyr::select(t2g, target_id, Geneid = ensembl_gene_id, Gene_Name = external_gene_name)

#Check the table "t2g" contains 'ens_gene'(Ensembl gene names) and target_id (associated transcripts from Ensembl).
head(t2g)

#reorder the columns
colnames(t2g) #Get column names

t2g <- t2g[, c(2, 3, 1)] # second column be the 1st, 3rd column be the 2nd and 1st column be the 3rd
head(t2g)


#merge two columns (counts and gene annotation table)

#counts1 <- merge(counts1,t2g, by="Geneid")
#head(counts1)
#dim(counts1)


# Import information about sample genotypes and conditions 
expdesign <- read.table("expdesign1.txt", header=T)

head(expdesign, n = 15)

dim(expdesign)


#check the experiment group we have

groups <- paste(expdesign$condition, expdesign$status, sep=".")

groups <- factor(groups) #group	vector containing the experimental group/condition for each sample(library)

table(groups) 

# 3. load required libraries
library("DESeq2")
library("edgeR")
library("limma")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("NMF")
library("EnhancedVolcano")
library("Rcpp")
library("pheatmap")

# let's get rid of some lowly expressed genes
#data_subset <- counts1[rowSums(counts1)>10,]
#dim(data_subset)



# Converting counts to DGEList object

class(counts1) #checking if is data frame

counts1 <- as.matrix(counts1) #convert df to matrix

y <- DGEList(counts1, samples=expdesign, group=expdesign$condition)

head(y)
#table(y)

#d <- DGEList(counts = data_subset, group=groups)
#head(d)


# 4. Filtering data
#using edgeR the raw counts are converted to CPM and log-CPM values using the cpm function.
#By default in edgeR, a gene should have CPM > 0.5 in at least 3 samples in opposite condition, the gene is removed

dim(y) #Checking before filtering step
#[1] 60664    13

keep <- rowSums(cpm(y) > 0.5) >= 2 #Identify genes with at least 0.5 cpm in at least 2 samples

table(keep)

y_filtered <- y[keep, , keep.lib.sizes = FALSE]


#data_subset1 <- cpm(counts1) #Calculate the Counts Per Million measure
#dim(data_subset1)


head(y_filtered)


dim(y_filtered) #Checking after filtering step
#[1] 31154    13
