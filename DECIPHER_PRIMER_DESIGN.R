rm(list = ls())

dir <- '/Users/cigom/metagenomics/db/fish_divergence/mafft/'
p <- '/Users/cigom/metagenomics/db/fish_divergence/mafft/'

setwd(dir)

rf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher.afa'
bf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher_vs_metafishgom_bold.afa'

nf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher_vs_metafishgom_ncbi.afa'

qr <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# specify the path to your sequence file:
fas <- paste0(dir, '/',bf)

# specify a path for where to write the sequence database
dbConn <- paste0(dir, '/DdbConn/')

library(DECIPHER)
library(tidyverse)
# OR create the sequence database in memory
dbConn <- dbConnect(SQLite(), ":memory:")

Seqs2DB(fas, "FASTA", dbConn, "")

# Defining groups ----
# Alternatively we could use the functions
# IdentifyByRank, FormGroups, or IdClusters to define groups.

ref <- readDNAStringSet(paste0(p, rf), format="fasta")
ncbi <- readDNAStringSet(paste0(p, nf), format="fasta")
bold <- readDNAStringSet(paste0(p, bf), format="fasta")

bold <- bold[!names(bold) %in% names(ref)]
ncbi <- ncbi[!names(ncbi) %in% names(ref)]

tax <- data.frame(name = names(bold))

tax %>% separate(sep = ' ', 
                 col = 1, 
                 into = c('id', 'Taxonomy')) %>%
  separate(sep = ';',col = 2, into = qr) %>%
  as_tibble() -> tax

target <- 'Albuliformes'
tax %>%
  filter(Order == target ) %>%
  select(id, Species) -> g

desc <- g[[2]]
  
names(bold) <- unlist(lapply(strsplit(names(bold), 
                                      " ", fixed=TRUE), 
                             function(x) return(x[1])))

seqgr <- bold[names(bold) %in% g[[1]],]
#DECIPHER::BrowseSeqs(seqgr)
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(seqgr, "XStringSet", dbConn, "")

tiles <- TileSeqs(dbConn, add2tbl="Tiles", minCoverage=1)
dim(tiles)

# Design primer
# Must install system("hybrid-min -V")
primers <- DesignPrimers(tiles, minCoverage=1, minGroupCoverage=1)

primers[1,]
names(primers)


temp_range <- 60:75

ps <- c("GCCTGGGGCTTTATTGGGAG", "AAAGCCCCAGGCTGACTTA") # forward and reverse
f <- function(temp) {
  CalculateEfficiencyPCR(ps, reverseComplement(DNAStringSet(ps)),
                         temp, P=4e-7, ions=.225)
}
efficiency <- matrix(unlist(lapply(temp_range, f)), ncol=2, byrow=TRUE)

plot(temp_range, efficiency[,1], ylim=c(0,1), 
       ylab="Hybridization Efficiency",
       xlab=expression(paste("Temperature (", degree, "C)", sep="")),
       type="l", lwd=2, col="Blue", main="Denaturation Plot")

lines(temp_range, efficiency[,2], col="Red", lwd=2)
abline(h=0.5, lty=2, lwd=2, col="Orange")
abline(v=64, lty=2, lwd=2, col="Green")

legend("topright", 
       legend=c("Forward Primer", "Reverse Primer", 
                "50% Efficiency", "Annealing Temperature"), 
       col=c("Blue", "Red", "Orange", "Green"),
         lwd=c(2, 2, 2, 2), lty=c(1, 1, 2, 2))

# Browser seqs
dna <- SearchDB(dbConn)
dbDisconnect(dbConn)
amplicon <- subseq(dna, 247, 348)
names(dna) <- desc

# only show unique sequences
u_amplicon <- unique(amplicon)

names(u_amplicon) <- names(amplicon)[match(u_amplicon, amplicon)]
amplicon <- u_amplicon
# move the target group to the top
w <- which(names(amplicon) == target)
amplicon <- c(amplicon[w], amplicon[-w])
BrowseSeqs(dna)


# test function!!!
minLength = 26 
maxLength = 27
maxTilePermutations = 10
minCoverage = 0.9


target <- SearchDB(dbConn, tblName = tblName, type = "DNAStringSet", 
                   identifier = identifier[k], processors = processors, 
                   verbose = FALSE, ...)

target <- SearchDB(dbConn)

tGaps <- TerminalChar(target)

a <- alphabetFrequency(target)
all(rowSums(a[, 1:15]) < maxLength)

consensus <- ConsensusSequence(target)

consensus <- ConsensusSequence(target)
pos <- which(strsplit(toString(consensus), "", fixed = TRUE)[[1]] != 
               "-")
l <- length(pos) - maxLength + 1

row_start <- 0
count <- 0
uw <- unique(width(target))
identifier = ""
k <- 0
tiles <- data.frame(row_names = (row_start + count + 
                                   1):(row_start + count + l * maxTilePermutations), start = I(integer(l * maxTilePermutations)), 
      end = I(integer(l * maxTilePermutations)), 
      start_aligned = I(integer(l * maxTilePermutations)), 
      end_aligned = I(integer(l * maxTilePermutations)), 
      misprime = I(logical(l * maxTilePermutations)), 
      width = I(rep(uw, l * maxTilePermutations)), 
      id = I(rep(identifier[k], l * maxTilePermutations)), 
      coverage = I(numeric(l * maxTilePermutations)), 
      groupCoverage = I(numeric(l * maxTilePermutations)), 
      target_site = I(character(l *maxTilePermutations)))
