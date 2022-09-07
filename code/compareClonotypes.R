rm(list = ls())
setwd("G:/My Drive/UCSD/research/Nathan E. Lewis/BCR/BCR_code_submission/DoubleCheckaFINAL")

#Load library
suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
library(magrittr)
library(reshape2)
library(dplyr)
library(ggalluvial)
#Read in 10X annotations
S1 <- read.csv("data/filtered_contig_annotations_NT-lsdB.csv")
S2 <- read.csv("data/filtered_contig_annotations_LAC-alum.csv")
S3 <- read.csv("data/filtered_contig_annotations_LAC-lsdB.csv")

S1$productive <- gsub('true','TRUE',S1$productive)
S2$productive <- gsub('true','TRUE',S2$productive)
S3$productive <- gsub('true','TRUE',S3$productive)


#Merge into a contig list
contig_list <- list(S1, S2, S3)

#Combine data together
combined <- combineBCR(contig_list, 
                       samples = c("NT","LAC","LAC"),
                       ID = c("lsdB","alum","lsdB2"))

###########################
#####compareClonotypes#####
###########################

#Remove list elements that contain all NA values
checkBlanks <- function(df, cloneCall) {
    for (i in seq_along(df)) {
        if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall]))) | 
            length(!is.na(df[[i]][,cloneCall])) == 0 | 
            nrow(df[[i]]) == 0) {
            df[[i]] <- NULL
        } else {
            next()
        }
    }
    return(df)
}

# This is to help sort the type of clonotype data to use
theCall <- function(x) {
    if (x == "gene") {
        x <- "CTgene"
    } else if(x == "nt") {
        x <- "CTnt"
    } else if (x == "aa") {
        x <- "CTaa"
    } else if (x == "gene+nt") {
        x <- "CTstrict"
    }
    return(x)
}

#compareClonotypes function
compareClonotypes <- function(df, cloneCall = "gene+nt", chain = "both", samples = NULL, 
                              clonotypes = NULL, numbers = NULL, graph = "alluvial",
                              exportTable = FALSE){
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    if (!is.null(numbers) & !is.null(clonotypes)) {
        stop("Make sure your inputs are either numbers or clonotype sequences.")
    }
    Con.df <- NULL
    for (i in seq_along(df)) {
        if (chain != "both") {
            df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
        }
        tbl <- as.data.frame(table(df[[i]][,cloneCall]))
        tbl[,2] <- tbl[,2]/sum(tbl[,2])
        colnames(tbl) <- c("Clonotypes", "Proportion")
        tbl$Sample <- names(df[i])
        Con.df <- rbind.data.frame(Con.df, tbl)
    }
    if (!is.null(samples)) {
        Con.df <- Con.df[Con.df$Sample %in% samples,] }
    if (!is.null(clonotypes)) {
        Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] }
    if (!is.null(numbers)) {
        top <- Con.df %>% top_n(n = numbers, wt = Proportion)
        Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] }
    if (nrow(Con.df) < length(unique(Con.df$Sample))) {
        stop("Reasses the filtering strategies here, there is not 
            enough clonotypes to examine.") }
    if (exportTable == TRUE) { return(Con.df)}
    
    plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                               stratum = Clonotypes, alluvium = Clonotypes, 
                               y = Proportion, label = Clonotypes)) +
        theme_classic() +
        theme(axis.title.x = element_blank())
    if (graph == "alluvial") {
        plot = plot +  geom_stratum() + geom_flow(stat = "alluvium")
    } else if (graph == "area") {
        plot = plot +
            geom_area(aes(group = Clonotypes), color = "black") }
    return(plot)
}

# remove NA
idx1 = which(is.na(combined$NT_lsdB$cdr3_aa2));   Label1 = paste0("Undefined",1:length(idx1))
idx2 = which(is.na(combined$LAC_alum$cdr3_aa2));  Label2 = paste0("Undefined",1:length(idx2))
idx3 = which(is.na(combined$LAC_lsdB2$cdr3_aa2)); Label3 = paste0("Undefined",1:length(idx3))
combined$NT_lsdB$CTaa  [idx1]= Label1
combined$LAC_alum$CTaa [idx2]= Label2
combined$LAC_lsdB2$CTaa[idx3]= Label3

# print figure 4A
num_clonotypes = 40
SampleSet = c("NT_lsdB","LAC_alum","LAC_lsdB2")

Results.DIR <- "results/"
if (!file.exists(Results.DIR)){dir.create(file.path(Results.DIR))}

fp = paste0("results/compareClonotypes_aa.lsdB2.N=",num_clonotypes,".pdf")
pdf(file=fp, width=15, height=8)
df.cast = list()
txt=paste0("clonotype=",num_clonotypes)
df <- compareClonotypes(combined, 
                        numbers = num_clonotypes, 
                        samples = SampleSet,
                        cloneCall="aa",
                        exportTable=TRUE)
df.cast[[as.character(num_clonotypes)]]<-dcast(df, Clonotypes~Sample, value.var = "Proportion")
print( compareClonotypes(combined, 
                         numbers = num_clonotypes, 
                         samples = SampleSet,
                         cloneCall="aa") + ggtitle(txt))
print(txt)
dev.off()

# write EXCEL file (Clonotypes & Sample) for figure 4A
fp2 = paste0("results/compareClonotypes_aa.lsdB2.N=",num_clonotypes,".xlsx")
write.xlsx(df.cast, file=fp2)

