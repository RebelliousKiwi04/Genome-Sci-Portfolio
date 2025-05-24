# 203311 Portfolio
### Chris Sandlant - 23003742

## Week 4 Portfolio
*Question is to determine whether GC content from either Oxford or Illumina reads is correlatd with quality scores*

For week 4 I have chosen to do the question with the Oxford Nanopore reads (Montana)

**Shell Script Code**


```console

# - g added to the options gets the GC content for the sequences
seqkit fx2tab -qlng montana-2021-29-09.fastq.gz > montana-fx2tabOut.txt

```

**R Script Code**

```R
##Read in initial data, and label columns for all of the new tables
montanaData <- read.table(file="montana-fx2tabOut.txt")

colnames(montanaData) <- c("sequence-name", "length", "gc content", "quality")


#Do initial correlation test - Week 4 Figure 1 below
cor.test(as.vector(montanaData$quality), as.vector(montanaData$`gc content`))

#Graphing in sorted line plot - Week 4 Figure 2 below
plot(sort(montanaData$'gc content'), sort(montanaData$quality), ty="l", col="purple", lwd=2, main="Montana GC Content vs Quality", ylab = "Sorted Montana Quality Data", xlab = "Sorted Montana GC Content")

```

### Week 4 Figure 1
![Week 4 Figure 1](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week4/Figure1.png)

The correlation test in this figure of -0.06 indicates a slightly negative relationship, but it is not nonlinear, as one increases, the other does tend to decrease slightly. This is reinforced by our sub 0.05 p-value, marking the correlation as statistically significant.

### Week 4 Figure 2
![Week 4 Figure 2](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week4/Figure2.png)

Our conclusions from figure 1 are further supported when analysing the graph, although it doesn't appear that way, this is because the seqkit quality scores are Phred quality scores, which are an ASCII quality score format, so therefore higher is not better, but worse. This indicates that the graph actually shows a negative relationship between GC content and quality for the montana sequence sample. Supporting our initial conclusion from figure 1, that as one increases, the other does tend to decrease.


## Week 5 Portfolio
*Question is to find a method to visualise read depth to determind if there is sufficient depth across the genome and note*

No shell script is required - depth text files were acquired from the lab itself

**R Script Code**

```R
#Read Depth information into R, and label columns 
kwazulu.depth <- read.table('kwazulu_depth.txt')
colnames(kwazulu.depth) <- c("Contig Name", "Contig Position", "Depth")

montana.depth <- read.table('montana_depth.txt')
colnames(montana.depth) <- c("Contig Name", "Contig Position", "Depth")

#Want to get mean depth - gives us a very early indication - Figure 1
mean(kwazulu.depth$Depth)
mean(montana.depth$Depth)


par(mfrow=c(1,2)) #Makes it possible to go dual panel

hist(montana.depth$Depth, col="powderblue", xlab="Montana Sequence Depth", main="Histogram of Montana sequence depth")

hist(kwazulu.depth$Depth, col="palegreen3", xlab ="Kwazulu-Natal Sequence Depth", main="Histogram of Kwazulu-Natal sequence depth") 

dev.off() #Switch off panel - resets effectively

#The figure generated from the code here is our figure 2 for Week 5

```

### Week 5 Figure 1
![Week 5 Figure 1](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week5/Figure1.png)

Description here

### Week 5 Figure 2
![Week 5 Figure 2](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week5/Figure2.png)

Description here


## Week 6 Portfolio
*Question is to use a new multifasta, align these new sequences, infer a phylogeny, assess the amount of certainty by performing bootstraps on the phylogeny, and plot the phylogeny with boostrapped values appearing on it*

Chose to do 1000 bootstrap cycles - as that was the number most often mentioned as the optimum for stable and accurate distribution

**Shell Script Code**

```console
#Get the multifasta
wget https://pjbiggs.github.io/Massey203311/Week6/data/covid_samples.fasta.gz

#Unzip it
gunzip covid_samples.fasta.gz

#Merge new file with our masked sequences
cat covid_samples.fasta montana-mask.fasta kwazulu-mask.fasta > all_genomes.fasta

#Perform the alignment
mafft --auto --reorder all_genomes.fasta >all_genomes.aln

#Now we generate the bootstrapped tree, -b command followed by number of bootstrap cycles
iqtree -b 1000 -s all_genomes.aln

#This generates our bootstrapped tree

```

**R Script Code**

```R
#From here we simply load the tree into R, and display it as we see fit
library(ape)

portfolio.tree <- read.tree(file="all_genomes.aln.treefile")

#Can be a bit squished - best results came fullscreen slightly zoomed out
plot.phylo(portfolio.tree, show.node.label = TRUE, tip.color="orange", node.color='darkblue', no.margin=TRUE)

#This generates our only figure for this week
```

### Week 6 Figure 1
![Week 6 Figure 1](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week6/Figure1.png)

Description here



## Week 8 Portfolio
*Week 8 has 2 parts to it, parts A and B*
> Part A involves the ``track`` dataframe, the task is to plot the dataframe in a way that shows whether the steps involved in rpocessing are adversely affecting one or more samples.

> Part B is to generate a sequence logo. Do this by generating a new fasta file of 50 sequences chosen at random from the reference fasta ``rdp_train_set_14.fa``, searching with ```msa`` and ``seqLogo`` to find a region of ~100-150bp in length that shows a high degree of sequence similarity and variability, then plot that region with the ``msa`` package showing the alignment and sequence logo.

There is no shell script requirement for this portfolio - All R based

**R Script Code**

```R
##PART A
#Copy track to new variable to avoid damaging original var
trackCopy <- track
#Proportionate data using prop.table
trackCopy <- prop.table(trackCopy)

#Now copy sample names to new column at the end and make row name numeric
trackCopy2 <- cbind(trackCopy, sample.names)

rownames(trackCopy2) <- 1:nrow(trackCopy2)

#Because i'm using a ggplot multi line plot - from R-Graph gallery - need to have a differently structured dataframe, so need to change how data is oriented, need to shrink it down to 3 columns
#effectively pivot it, only way i could think to do this was this way - really really messy, alternative was a library I found - declared in AI declaration - I chose to use this method instead
outputFrame <- data.frame() #Make new dataframe to fill
for (row in 1:nrow(as.matrix(trackCopy2))) {#Iterate through rows
    #print(trackCopy2[row,7])
    appendSample <- trackCopy2[row,7]#Grab sample name
    for(col in 1:6) {#Iterate through only the first 6 columns (we dont want to grab the sample names)
      appendStage <-colnames(trackCopy2)[col] #Grab stage name
      #print(trackCopy2[row,col])
      appendValue <-trackCopy2[row,col] #Grab value
      #print(colnames(trackCopy2)[col])
      appendList <- list(appendSample, appendStage, appendValue) #Create list
      outputFrame <- rbind(outputFrame, appendList)#Bind list to dataframe as new row
      
    }
}
#Rename columns again
colnames(outputFrame) <- c("sample.names", "stage", "value")

#Do ggplot line plot, with stage on the x axis, and vlaues on y axis, coloured by samples
ggplot(outputFrame, aes(x = stage, y = value, group = sample.names, color = sample.names)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +
  theme_light() +
  labs(
    title = "Proportionated Track Analysis",
    x = "Processing Stage",
    y = "Read Count",
    color = "Sample Name"
  )#Figure 1


##PART B
library(FastaUtils)
library(msa)
library(seqLogo)

#This gets us our random 50 sequences to a  new fasta file - sub50_rdpset.fa
fasta.sample(infile="rdp_train_set_14.fa", nseq=50, file.out="sub50_rdpset.fa")

#Create variable from subset
seqs16S <- readDNAStringSet("sub50_rdpset.fa")

#Run sequence alignment
my16SAlignment <- msa(seqs16S, "ClustalOmega")

#Print sequence to pdf for region selection
msaPrettyPrint(my16SAlignment, output="pdf", showNames="none", file = "ourSet.pdf",
               showLogo="top", askForOverwrite=FALSE, verbose=FALSE)
               #Figure 2 comes from here

#Identified region of high similarity and variability is 1000-1150
#This creates new object  using only positions 1000-1150
smallAlign <- DNAStringSet(my16SAlignment, start=1000, end=1150)
#Creates consensus matrix usingn seqLogo
posWeightMat <- consensusMatrix(smallAlign, baseOnly=TRUE)
#Remove other row
rowname_to_remove <- ("other")
posWeightMat2 <- posWeightMat[!(row.names(posWeightMat) %in% rowname_to_remove),]
#proportionate and make seqLogo
posWeightMat2 <- prop.table(posWeightMat2, 2)
p <- makePWM(posWeightMat2)
seqLogo(p, xfontsize = 8)#SeqLogo is Figure 3

```
### Week 8 Figure 1
![Week 8 Figure 1](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week8/Figure1.png)

Description here

### Week 8 Figure 2
![Week 8 Figure 2](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week8/Figure2.png)

Description here

### Week 8 Figure 3
![Week 8 Figure 3](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week8/Figure3.png)

Description here


## Week 9 Portfolio
*Week 9 has 2 parts A, and B*
Part A - Chose Combination 2 - ``500k_Cutoff1.txt`` and taxonomic classification ``__Rhodobacter``
Part B - Chose Combination 2 - ``500k_Cutoff1.txt``

> Part A involves trimming out the constant part of the taxonomic name to reduce the length of classification name and visualising the data in a heatmap

> Part B Involves visualising the chosen file in an UpSet plot, to show how the number of times a taxonomic result varies bbased on searching algorithm and database

There is no shell script requirement for this portfolio - All R based

**R Script Code**

```R
library(UpSetR)
library(pheatmap)
library(data.table)

#Firstly we will do Part A

#Read in the tab delimited data file
cutDataFile <- read.delim("500k_Cutoff1.txt", header=TRUE, sep="\t")
#Remove columns
cutDataFileRemovedCols <- subset(cutDataFile, select=-c(averVal, COV))
cutDataFileRemovedCols_noUC <- cutDataFileRemovedCols[-c(1), ]

rownames(cutDataFileRemovedCols_noUC) <- cutDataFileRemovedCols_noUC[,7]
cutDataFileRemovedCols_noUC2 <- cutDataFileRemovedCols_noUC[,-7]
#Show only taxa that contain Rhodobacter
onlyTaxaOfInterest <- subset(cutDataFileRemovedCols_noUC2, rownames(cutDataFileRemovedCols_noUC2) %like% "__Rhodobacter__")

##this line gsubs out the prefix
rownames(onlyTaxaOfInterest) <- gsub("root__cellular__organisms__Bacteria__Proteobacteria__Alphaproteobacteria__Rhodobacterales__Rhodobacteraceae__Rhodobacter__", "", rownames(onlyTaxaOfInterest))

#Generate heatmap - Figure 1
pheatmap(onlyTaxaOfInterest)


#Now onto Part B

#Removes unwanted columns from our dataset
cutDataFileUpSet <- subset(cutDataFile, select=-c(averVal, COV, taxonomy))
#Removes non-0 values to 1
cutDataFileUpSet[cutDataFileUpSet != 0] <- 1
#Add our taxonomy back into the set
cutDataFileUpSetTaxa <- cbind(cutDataFileUpSet, cutDataFile$taxonomy)
colnames(cutDataFileUpSetTaxa)[7] <- "taxonomy"

#Draw plot - Figure 2
setOrder = c("Greedy_Nr", "Greedy_NrEuk", "Greedy_RefOnly", "MEM_Nr", "MEM_NrEuk", "MEM_RefOnly")
upset(cutDataFileUpSetTaxa, order.by = "freq", decreasing = TRUE,sets = setOrder, keep.order = TRUE)

```

### Week 9 Figure 1
![Week 9 Figure 1](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week9/Figure1.png)

Description

### Week 9 Figure 2
![Week 9 Figure 2](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week9/Figure2.png)

Description

## Week 10 Portfoliio
*Week 10 again has 2 parts, working with the portfolioPhy phyloseq object we made at the end of the lab, making a new phyloseq object subsetted at the taxonomic level of kingdom, containing only Bacteria*
> Part A involves plotting alpha diversity values with ``Sample`` on the x axis, colouring by ``env_broad_scale`` with the ``Shannon`` diversity measure. And a second plot, plotting the values for ``env_local_scale`` on the x-axis colouring by ``env_broad_scale`` with both the ``Shannon`` and ``Simpson`` diversity measures.

> Part B involves plotting an single barplot, according to the following
- Make a new object of the top 100 taxa
- Transform the sample conts to a fraction incorporating the code ``function(OTU OTU/sum(OTU))
- Prune the taxa
- Plot the barplot by ``env_local_scale`` on the x-axis with colouring by ``Genus``~and ``facet_wrap`` on ``env_broad_scale``

There is no shell script for this week's portfolio - all R based

**R Script code**
```R
#Create new subsetted phyloseq object
portfolio.bac <- subset_taxa(portfolioPhy, Kingdom=="Bacteria")

#Part A First graph - Figure 1
plot_richness(portfolio.bac, x="Sample", measures=(c("Shannon")), color="env_broad_scale")

#Part A second graph - Figure 2
plot_richness(portfolio.bac, x="env_local_scale", measures=c("Shannon", "Simpson"), color="env_broad_scale")

#Part B
#Get top 100 taxa in descending order
portfolioTop100 <- names(sort(taxa_sums(portfolio.bac), decreasing=TRUE))[1:100]
#Transform by function
portfolio.top100 <-transform_sample_counts(portfolio.bac, function(OTU) OTU/sum(OTU))
#Prune
portfolio.top100 <-prune_taxa(portfolioTop100, portfolio.top100)
#Plot barplot - Figure 3
plot_bar(portfolio.top100, x="env_local_scale", fill="Genus") +facet_wrap(~env_broad_scale, scales="free_x")

```

### Week 10 Figure 1
![Week 10 Figure 1](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week10/Figure1.png)

Descrition

### Week 10 Figure 1
![Week 10 Figure 2](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week10/Figure2.png)

Description

### Week 10 Figure 1
![Week 10 Figure 3](https://rebelliouskiwi04.github.io/Genome-Sci-Portfolio/Figures/Week10/Figure3.png)

Description


## Week 11 Portfolio
*Week 11s portfolio requires the calculation of sequence depth for 6 bam files, that then need to be plotted in such a way that the distributions are easy to compare, and thus easy to check for problematic samples*

**Shell script**

```console

#Shell script side is simple, just getting depth information for the sequences

#UHR Samples
samtools depth UHR_Rep1.sort.bam > UHR_Rep1.depth.txt

samtools depth UHR_Rep2.sort.bam > UHR_Rep2.depth.txt

samtools depth UHR_Rep3.sort.bam > UHR_Rep3.depth.txt

#HBR Samples
samtools depth HBR_Rep1.sort.bam > HBR_Rep1.depth.txt

samtools depth HBR_Rep2.sort.bam > HBR_Rep2.depth.txt

samtools depth HBR_Rep3.sort.bam > HBR_Rep3.depth.txt

```

**R Script code**

```R

#First load depth information in as table
UHR_Rep1 <- read.table(file="UHR_Rep1.depth.txt")
UHR_Rep2 <- read.table(file="UHR_Rep2.depth.txt")
UHR_Rep3 <- read.table(file="UHR_Rep3.depth.txt")

HBR_Rep1 <- read.table(file="HBR_Rep1.depth.txt")
HBR_Rep2 <- read.table(file="HBR_Rep2.depth.txt")
HBR_Rep3 <- read.table(file="HBR_Rep3.depth.txt")

#Name columns - reference Week 5 - same layout
colnames(UHR_Rep1) <- c("Contig Name", "Contig Position", "Depth")
colnames(UHR_Rep2) <- c("Contig Name", "Contig Position", "Depth")
colnames(UHR_Rep3) <- c("Contig Name", "Contig Position", "Depth")

colnames(HBR_Rep1) <- c("Contig Name", "Contig Position", "Depth")
colnames(HBR_Rep2) <- c("Contig Name", "Contig Position", "Depth")
colnames(HBR_Rep3) <- c("Contig Name", "Contig Position", "Depth")

#Just need to decide on final plotting format
#Either hist,density, or violin
#All have their problems, could go with all 3



```