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

