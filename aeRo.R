####INSTALLATION####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")


####LIBRARY####
library(dada2); packageVersion("dada2")

####


####SET WORKING DIRECTORY####
path <- "~/Desktop/aeroinfo/fastq" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
setwd(path)


####GET LISTS####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



####INSPECT QUALITY####
plotQualityProfile(fnFs[1:2]) #outputs a quality plot
plotQualityProfile(fnRs[1:2]) #outputs a quality plot


####FILTERING####
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


####LEARN ERROR RATES####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) #error plot Forward
plotErrors(errR, nominalQ=TRUE) #error plot Reverse


####SAMPLE INFERENCE####
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #forward algorithm
dadaRs <- dada(filtRs, err=errR, multithread=TRUE) #reverse algorithm

#inspect the dada class object
dadaFs[[1]]


####MERGE PAIRED END READS####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


####CONSTRUCT SEQUENCE TABLE####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


####REMOVE CHIMERAS####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab) #inspect chimeric ratio


####TRACKING READS THROUGH PIPELINE####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


####ASSIGN TAXONOMY####
#what database should I use for this? Check out c20a paper and find out
#greengenes/silva?
#going to go with silva database
#genus level ID:
taxa <- assignTaxonomy(seqtab.nochim, "/Users/Tristan/Desktop/aeroinfo/fastq/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#species level ID:
taxa <- addSpecies(taxa, "/Users/Tristan/Desktop/aeroinfo/fastq/silva_species_assignment_v132.fa.gz")
#inspect
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) #looks good


####HANDOFF TO PHYLOSEQ####
#Libraries
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

#####IMPORT METADATA####
#mdimport is metadata that is 'messy' because rownames are numbers, not the names of your samples! Need to clean them up
mdimport <- read.table(file = "/Users/Tristan/Desktop/aeroinfo/fastq/GLDS-170_sample_metadata.csv", header=TRUE, sep = ",")
mdimport <- as.data.frame(md)
md <- data.frame(mdimport[,-1], row.names=mdimport[,1]) #cleaves off the indices that R mistakenly reads as rownames
md #md is the cleaned up metadata
identical(rownames(seqtab.nochim), rownames(md)) #sanity check to ensure that rownames match

####PHYLOSEQ OBJECT####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(md), 
               tax_table(taxa))
#Save as CSV
write.csv(seqtab.nochim, file = "ps.csv")

####Abbreviating ASVs; Creating refseq####
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

####DECONTAM####
#install decontam
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
head(sample_data(ps)) #not being used
####Inspect Library Sizes####
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Characteristics..Sample.Category)) + geom_point()


####Identify Contaminants - Prevalance####
sample_data(ps)$is.neg <- sample_data(ps)$Characteristics..Sample.Category == "Ground"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

