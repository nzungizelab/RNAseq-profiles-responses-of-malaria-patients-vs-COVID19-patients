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


#among various normalisation methods (TMM, BH, RLE, upperquartile UQ) we will use  TMM 

#These normalization are performed using either DESeq (RLE) or edgeR (TC, RPKM, UQ, TMM).

#TMM: estimating relative RNA production levels from RNA-seq data
y_filtered <- calcNormFactors(y_filtered, method="TMM")  #Calculate Normalization Factors to Align Columns of a Count Matrix 

head(y_filtered$samples, n=15)


#design matrix
design <- model.matrix(~0+groups)


colnames(design) <- levels(groups)

design


y_filtered <- estimateDisp(y_filtered, design, robust=TRUE)

y_filtered$common.dispersion
#[1] 0.5400794

#edgeR uses the Cox-Reid profile-adjusted likelihood method in estimating dispersion
y_filtered <- estimateCommonDisp(y_filtered,verbose=TRUE)
#Disp = 0.51809 , BCV = 0.7198 



fit <- glmQLFit(y_filtered, design, robust=TRUE) #estimated values of the GLM coefficients for each gene.

head(fit$coefficients)

summary(fit$df.prior)


COVID19vsHealthy <- makeContrasts( COVID19.patient-Healthy.donor, levels=design)

MalariavsHealthy <- makeContrasts( Malaria.patient-Healthy.donor, levels=design)

COVID19vsMalaria <- makeContrasts( COVID19.patient-Malaria.patient, levels=design)


#DE in each group using the QL F-test.

#DE in Covid19 patient
COVID19vsHealthy_Tab1 <- glmQLFTest(fit, contrast=COVID19vsHealthy, coef = 2)


topTags(COVID19vsHealthy_Tab1, n= 5) # table of top 5 DEG

COVID19vsHealthy_top.gene <- topTags(COVID19vsHealthy_Tab1)  # table of top DEG
COVID19vsHealthy_top.gene




COVID19vsHealthy_de <- decideTestsDGE(COVID19vsHealthy_Tab1, p.value=0.01) # pvalue using Wald test method


summary(COVID19vsHealthy_de) 
# Down o gene and Up 0 gene

# adjust p-values and assign the result to our table
#we consider a fraction of 10% false positives acceptable, therefore all genes with an adjusted p value below 10% = 0.1 as significant. 
COVID19vsHealthy_Tab1$table$padj <- p.adjust(COVID19vsHealthy_Tab1$table$PValue, method="BH")

sum(COVID19vsHealthy_Tab1$table$padj < 0.1)
# 0 gene

head(COVID19vsHealthy_Tab1$table$padj)

topTags(COVID19vsHealthy_Tab1)



# highlight DE genes in orange
points(COVID19vsHealthy_Tab1$table$logCPM[COVID19vsHealthy_de], COVID19vsHealthy_Tab1$table$logFC[COVID19vsHealthy_de], col="orange")



# How many genes look significant?
sum(COVID19vsHealthy_Tab1$table$PValue < 0.01)
# 2866 genes

sum(COVID19vsHealthy_Tab1$table$padj < 0.1)
# padj < 0.01 = 0 genes;  padj < 0.1 = 2495


# How many genes show 2-fold enrichment?
sum(COVID19vsHealthy_Tab1$table$PValue < 0.01 & COVID19vsHealthy_Tab1$table$logFC > 1)
#251 genes

sum(COVID19vsHealthy_Tab1$table$padj < 0.1 & COVID19vsHealthy_Tab1$table$logFC > 1)
# # padj < 0.01 = 0 genes;  padj < 0.1 = 188


#DE in Malaria infection 
MalariavsHealthy_Tab2 <- glmQLFTest(fit, contrast=MalariavsHealthy, coef = 2)

topTags(MalariavsHealthy_Tab2, n= 5) # table of top 5 DEG

MalariavsHealthy_top.gene <- topTags(MalariavsHealthy_Tab2)  # table of top DEG
MalariavsHealthy_top.gene

MalariavsHealthy_de <- decideTestsDGE(MalariavsHealthy_Tab2, p.value=0.01) 


summary(MalariavsHealthy_de) 
#Down 172 genes and Up 9 genes

# adjust p-values and assign the result to our table
MalariavsHealthy_Tab2$table$padj <- p.adjust(MalariavsHealthy_Tab2$table$PValue, method="BH")

sum(MalariavsHealthy_Tab2$table$padj < 0.1)
#padj < 0.01 = 181 genes;  padj < 0.1 = 1084 genes

head(MalariavsHealthy_Tab2$table$padj)

topTags(MalariavsHealthy_Tab2)

## Heatmap for Tab3_overlap

library(gplots)
library(RColorBrewer)

#Reading in data and transform it into matrix format
data <- read.csv("cov_mal.de_logFC.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
cov_mal_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-3 into a matrix
rownames(cov_mal_data) <- rnames                  # assign row names

#Customizing and plotting the heat map

# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("red", "yellow", "darkgreen"))(n = 62) #"red","yellow","green"

my_palette <-colorRampPalette(brewer.pal(3,"RdBu"))(61)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

# creates a 5 x 5 inch image
png(paste0(outpath,"heatmap.png"), # create and save PNG for the heat map 
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

# Heatmap performs reordering using clusterization
heatmap.2(cov_mal_data,
          #cellnote = cov_mal_data,  # same data set for cell labels
          #main = "Overlapped genes expressed", # heat map title
          #xlab="samples with infections", 
          #ylab="Overlapped genes expressed", 
          notecol = "black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(7,10),     # widens margins around plot #c(12,10)
          lwid = c(0.2,4),       # c(0.2,5)
          lhei = c(0.2,8),         #c(0.2,5)
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram= "row",     # only draw a row dendrogram
          #dendrogram="both",
          #Colv = FALSE, 
          #scale = "none",
          #scale = "row",  # #scale by row
          scale = "column",
          #key.xlab = "Abundance",
          Rowv = TRUE,  #TRUE or NA
          #cexRow = 2,
          #cexCol = 2,
          Colv= "NA"  # turn off column clustering
          )            

dev.off()               # close the PNG device

