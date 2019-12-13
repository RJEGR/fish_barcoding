
rm(list = ls())

# =======
# Library
# =======
.bioc_packages <- c("Biostrings", "ShortRead") # "fastqcr", ""
.inst <- .bioc_packages %in% installed.packages()

# if(any(!.inst)) {
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install(.bioc_packages[!.inst], ask = F)
# }

# Load packages into session, and print package version

sapply(.bioc_packages, require, character.only = TRUE)

# Functions:
# ==========

coordAlign <- function(Pattern,subject) {
  # Find merge region size
  # Make seed aligment between paired reads
  # First use a fixed substitution matrix
  
  class(Pattern)
  class(subject)
  
  mat <- nucleotideSubstitutionMatrix(
    match = 1, mismatch = -1, baseOnly = TRUE)
  
  pwa <- pairwiseAlignment(Pattern, subject, 
                           type = "local", 
                           substitutionMatrix = mat,
                           gapOpening = 0, gapExtension = 1,
                           scoreOnly = FALSE)
  
  alignment <- c(alignedPattern(pwa))
  
  # names(alignment) == names(Pattern)
  
  # alignment <- as.character(alignment)
  # 
  # matchPattern <- matchPattern(alignment, 
  #                              as.character(Pattern))
  # 
  # 
  # out <- vector("numeric", length = 4)
  # 
  #       
  # out[1] <- width(Pattern)
  # out[2] <- start(matchPattern)
  # out[3] <- end(matchPattern)
  # out[4] <- width(matchPattern)
  #out[5] <- names(Pattern)
    
  # length = width(Pattern),
    # start = start(matchPattern),
    # end = end(matchPattern),
    # width = width(matchPattern),
    # names = names(Pattern))
  
  return(alignment)
  
}

pattStats <- function(x) {
  x <- vector("numeric", length = 4)
  
  x[1] <- width(Pattern)
  x[2] <- start(matchPattern)
  x[3] <- end(matchPattern)
  x[4] <- width(matchPattern)
  return(x)
  
}

# length = width(Pattern),
# start = start(matchPattern),
# end = end(matchPattern),
# width = width(matchPattern),
# names = names(Pattern))
# =====
# Files
# =====
# Quality
qc1 <- "/Users/cigom/metagenomics/MG_18S/run04_20170418_18S/fastqc/04-X04-A10-18S-AMB_S7_L001_R1_001_fastqc.zip"
qc2 <- "/Users/cigom/metagenomics/MG_18S/run04_20170418_18S/fastqc/04-X04-A10-18S-AMB_S7_L001_R2_001_fastqc.zip"
# Reads
fnFs <- "/Users/cigom/metagenomics/MG_18S/run04_20170418_18S/04-X04-A8-18S-AMB_S6_L001_R1_001.fastq.gz" 

fnRs <- "/Users/cigom/metagenomics/MG_18S/run04_20170418_18S/04-X04-A8-18S-AMB_S6_L001_R2_001.fastq.gz" 

source('https://raw.githubusercontent.com/RJEGR/metagenomics/master/plotQualityProfile.R')

# plotQP(c(fnRs, fnFs))

# Find merge region size
# Make seed aligment between paired reads
# First load dataset

s1 <- readDNAStringSet(fnFs, format="fastq")
s2 <- readDNAStringSet(fnRs, format="fastq")
s2 <- Biostrings::reverseComplement(s2)

# f <- FastqStreamer(fnFs, n = 1e+06) # sampling records from fastq
seq1 <- s1[1:100]
seq2 <- s2[1:100]

coordAlign(seq1[1], seq2[1])

f1 <-mapply(coordAlign, seq1, seq2)

# Prepare some containers with correct length for 
# storing the result and change print() to assignment statement

f1 <- list()
#f2 <- list()

for (i in 1:length(seq1)) { f1[[i]] <- coordAlign(seq1[i], seq2[i])}
#for (i in 1:length(seq2)) { f2[[i]] <- coordAlign(seq2[i], seq1[i])}

dim(f1 <- data.frame(do.call(rbind, f1)))
#dim(f2 <- data.frame(do.call(rbind, f2)))

names <- c("length", "start", "end", "width")
names(f1) <- names
#names(f2) <- names

# cualquier alinemiento, seq1 vs seq2 o seq2 vs seq1 te regresa el tamano del empalme correspondiente para esas secuencias pareadas. 

# Existe alguna relacion entre el tamano del empalme vs el tamano de calidad en la region del empalme?

hist(f1$width)
coordAlign(seq1[10], seq2[10])

#seed1 <- subseq(s1, 1,30)[1]
#seed2 <- subseq(s2, 1,30)[1]

mat <- nucleotideSubstitutionMatrix(
  match = 1, mismatch = -1, baseOnly = TRUE)

pwa <- pairwiseAlignment(seq1, seq2, type = "local", 
                  substitutionMatrix = mat,
                  gapOpening = 0, gapExtension = 1,
                  scoreOnly = FALSE)


aligned_f1 <- c(alignedPattern(pwa))
#aligned_f2 <- c(alignedSubject(pwa))

aligned_f1 <- as.character(aligned_f1)
# aligned_f2 <- as.character(aligned_f2)

f1 <- list()

for (i in 1:length(aligned_f1)) { f1[[i]] <- matchPattern(aligned_f1[i], 
                                                          as.character(seq1)[i]) }

f1 <- data.frame(do.call(rbind, f1))


barplot(table(width(f1)))


# Average read quality vs read length

fastq <- data.frame(qc_read(qc2, 
                            modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)


# ================
# Using shortReads
# ================
# https://kasperdanielhansen.github.io/genbioconductor/html/ShortRead.html
library(ShortRead)

reads <- readFastq(fnFs)

# rl <- width(reads)
# qq <- as(quality(reads), "matrix")
# qq <- apply(qq, 1, mean)
# 
# dt <- data.frame(rl, qq)
# 
# d <- ggplot(dt, aes(rl, qq))
# d + geom_hex(bins = 30)





srqa <- qa(fnFs, n = 5e+05)

df <- srqa[["perCycle"]]$quality

d <- ggplot(df, aes(Cycle, Score))
d + geom_hex(bins = 30)

rc <- sum(srqa[["readCounts"]]$read)

means <- rowsum(df$Score * df$Count, df$Cycle)/ rowsum(df$Count, df$Cycle)
# quality

quality <- quality(reads)

encoding(quality)

# convert qualities into standard 0-40 integer quality scores
qquality <- as(quality(reads), "matrix")

hist(rowMeans(qquality[sample(nrow(qquality), 25000),],na.rm=TRUE), 
     main="", xlab="")

boxplot(qquality[sample(nrow(qquality), 25000), 1:100],out.cex=0.5)

mung <- function(i) { 
  lower <- i  
  upper <- i + 4  
  if (upper > ncol(qquality)) {    
    upper <- ncol(qquality) } 
  return(rowMeans(qquality[, seq(lower, upper)], na.rm=TRUE))
  }

sequence <- seq(1, ncol(qquality), 5)
groupedQuality <- as.data.frame(sapply(sequence, mung))

colnames(groupedQuality) <- as.character(sequence)

boxplot(groupedQuality[sample(nrow(groupedQuality), 25000),], outcex=0.5)

# Sampling
sampleSize <- 1:100
qquality <- qquality[sampleSize,]

plot(apply(qquality, 1, mean), type = 'l')

# [1] "Base"             "Mean"             "Median"          
# [4] "Lower.Quartile"   "Upper.Quartile"   "X10th.Percentile"
# [7] "X90th.Percentile"



summary()
as(quality(reads), "matrix")[1,1:151]
plot(as(quality(reads), "matrix")[1,1:151], type = 'l')
