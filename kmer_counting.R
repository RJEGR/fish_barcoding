rm(list = ls())

LCO1490 <- DNAString("GGTCAACAAATCATAAAGATATTGG")
HCO2198 <-  DNAString("TAAACTTCAGGGTGACCAAAAAATCA")

rf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher.afa'
bf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher_vs_metafishgom_bold.afa'
nf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher_vs_metafishgom_ncbi.afa'

target <- 'metafishgom_DB_target.tsv'

p <- '/Users/cigom/metagenomics/db/fish_divergence/mafft/'
wd <- '/Users/cigom/metagenomics/db/fish_divergence/'
setwd(wd)
# Load files

ref <- readDNAStringSet(paste0(p, rf), format="fasta")
ncbi <- readDNAStringSet(paste0(p, nf), format="fasta")
bold <- readDNAStringSet(paste0(p, bf), format="fasta")


# Remove the reference from the aligment

bold <- bold[!names(bold) %in% names(ref)]
ncbi <- ncbi[!names(ncbi) %in% names(ref)]

# Simulate Amplification of DNA by PCR
primers <- c(as.character(LCO1490), as.character(HCO2198))
myDNAStringSet <- DECIPHER::RemoveGaps(bold)

y <- AmplifyDNA(primers,
           myDNAStringSet,
           maxProductSize = 650,
           annealingTemp = 55,
           P = 4e-7,
           ions = 0.2,
           includePrimers = TRUE,
           minEfficiency = 0.001)


OrientNucleotides()
TrimDNA()

















# names(gom)[names(gom) %in% morf] <- morf_names
# names(gom)[names(gom) %in% wr[3:6]] <- rank[3:6]



#kmer Fast K-mer Counting and Clustering for Biological Sequence Analysis.

# kmer R package
# https://github.com/shaunpwilkinson/kmer
# https://cran.r-project.org/web/packages/kmer/kmer.pdf
# devtools::install_github("shaunpwilkinson/kmer") 

rm(list = ls())


library("kmer")
data(woodmouse, package = "ape")
ape::as.character.DNAbin(woodmouse[1:5, 1:5])

woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]


### Compute the full distance matrix and print the first few rows and columns

woodmouse.kdist <- kdistance(woodmouse, k = 6)
print(as.matrix(woodmouse.kdist)[1:7, 1:7], digits = 2)


### Compute and print the embedded distance matrix
suppressWarnings(RNGversion("3.5.0"))
set.seed(999)
seeds <- sample(1:15, size = 3)
# Convert sequences to vectors of distances to a subset of seed sequences.

woodmouse.mbed <- mbed(woodmouse, seeds = seeds, k = 6)

print(woodmouse.mbed[,], digits = 2)

## compute pairwise distance matrices
dist1 <- ape::dist.dna(woodmouse, model = "K80") 
dist2 <- kdistance(woodmouse, k = 7) 

## build neighbor-joining trees
phy1 <- ape::nj(dist1)
phy2 <- ape::nj(dist2)

## rearrange trees in ladderized fashion
phy1 <- ape::ladderize(phy1)
phy2 <- ape::ladderize(phy2)

## convert phylo objects to dendrograms
dnd1 <- as.dendrogram(phy1)
dnd2 <- as.dendrogram(phy2)

## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 5)

# 
# Installation of gkmSVM package:

#BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19.masked'))
#install.packages('ROCR')

# git clone https://github.com/mghandi/gkmSVM.git
# R CMD INSTALL gkmSVM
#Or  install.packages('gkmSVM')

library(gkmSVM)

#
