# IF not already installed, install the following packages. IF already installed, skip to next chunk
```{r}
##install.packages("biostrings")
##install.packages("dada2")
##install.packages("phyloseq")
##install.packages("ggplot2")
```
#load required libraries
```{r}
library(biostrings)
library(dada2)
library(phyloseq)
library(ggplot2)
```
#MANUAL CHANGE REQUIRED. make an object of working directory - so don’t have to type out the path again
```{r}
##path <-"CHANGEME_DIRECTORYPATH"
```
#list files in the "path" directory to doublecheck
```{r}
list.files(path)
```
#make a fnFs object of forward reads and a fnRs file of reverse reads
#MANUAL CHANGE REQUIRED for file pattern
```{r}
##fnFs <- sort(list.files(path, pattern="CHANGEME", full.names =TRUE))
##fnRs <- sort(list.files(path, pattern="CHANGEME", full.names =TRUE))
```
#check the quality of the data - look at green (mean quality score) line
```{r}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
```{r}
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```
#filter and trim the data. Place filtered files in /filtered/ subdirectory and create objects
```{r}
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```
#MANUAL CHANGE REQUIRED X Y. truncLen=c(x,y) cut forward reads at x basepairs and cut reverse reads at y basepairs according to quality scores
```{r}
##out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(x,y),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
head(out)
```
#visualize and learn the error rate
```{r}
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
```
#visualize sample interference
```{r}
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
```
#merge paired reads
```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```
# Inspect the merger data.frame from the first sample
```{r}
head(mergers[[1]])
```
#construct a sequence table seqtab
```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```
# Inspect distribution of sequence lengths
```{r}
table(nchar(getSequences(seqtab)))
```
# remove chimeras and make an object seqtab.nochim
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```
#look at the number of reads that made it through each step in the pipeline
```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
```
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
```{r}
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```
#first download and place taxonomy learning files into the working directory
#MANUAL CHANGE REQUIRED. assign taxanomy. irectory
```{r}
##taxa <- assignTaxonomy(seqtab.nochim,"CHANGEmeDIRECTORYPATH_LEARNINGFILE", multithread=FALSE)
```
#look into the taxa
```{r}
taxa.print <- taxa 
```
#Removing sequence rownames for display only
```{r}
rownames(taxa.print) <- NULL
head(taxa.print)
```
#MANUAL CHANGE REQUIRED. download the data into a csv file
```{r}
##write.csv(taxa, file = "CHANGEME/taxa.csv")
##write.csv(seqtab.nochim, file="CHANGEME/seqtab_nochim.csv")
```
#import data into phyloseq
```{r}
seqtab.nochim<- as.matrix(seqtab.nochim)
taxa<- as.matrix(taxa)
```
#take rownames as an object - turn into dataframe - put rownames as rownames 
```{r}
samples.out<- rownames(seqtab.nochim)
samdf<- data.frame(samples.out)
rownames(samdf)<- samples.out
```
#import data (abundance, taxa, sample names) into phyloseq as an object - ps, otu_table = abundance table in phyloseq
```{r}
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
```
#give ASVs simpler names i.e. ASV1, ASV2...
```{r}
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```
#plot with ggplot: need to take data out of phyloseq format and into a dataframe using psmelt function
```{r}
ps.table<- psmelt(ps)
```
#make the phylum column into a factor so that ggplot can process it
```{r}
ps.table$Phylum<- factor(ps.table$Phylum)
```
#make a ggplot now for phylum relative abundance - for bar plots
```{r}
ggplot(data = ps.table, mapping = aes(x=Sample,y=Abundance)) + geom_bar(aes(fill=Phylum), stat="identity", position = "fill")+labs(x="Sample", y="Abundance")
```
