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

