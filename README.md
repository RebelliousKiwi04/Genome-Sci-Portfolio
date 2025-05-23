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
![Week 4 Figure 1](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week4/Figure1.png)

Description text here

### Week 4 Figure 2
![Week 4 Figure 2](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week4/Figure2.png)

Description text here

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
![Week 5 Figure 1](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week5/Figure1.png)

Description here

### Week 5 Figure 2
![Week 5 Figure 2](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week5/Figure2.png)

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
![Week 6 Figure 1](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week6/Figure1.png)

Description here



## Week 8 Portfolio


IN PROGRESS - kinda a disaster at the moment


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
![Week 9 Figure 1](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week9/Figure1.png)

Description

### Week 9 Figure 2
![Week 9 Figure 2](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week9/Figure2.png)

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
![Week 10 Figure 1](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week10/Figure1.png)

Descrition

### Week 10 Figure 1
![Week 10 Figure 2](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week10/Figure2.png)

Description

### Week 10 Figure 1
![Week 10 Figure 3](https://github.com/RebelliousKiwi04/Genome-Sci-Portfolio/blob/master/Figures/Week10/Figure3.png)

Description


## Week 11 Portfolio

IN PROGRESS