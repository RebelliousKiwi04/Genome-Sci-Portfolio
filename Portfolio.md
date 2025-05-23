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
![Week 4 Figure 1](Figures\Week4\Figure1.png)

Description text here

### Week 4 Figure 2
![Week 4 Figure 2](Figures\Week4\Figure2.png)

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

## Week 5 Figure 1
![Week 5 Figure 1](Figures\Week5\Figure1.png)

Description here

## Week 5 Figure 2
![Week 5 Figure 2](Figures\Week5\Figure2.png)

Description here



