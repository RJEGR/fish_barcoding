rm(list = ls())

options(stringsAsFactors = FALSE)

.bioc_packages <- c("Biostrings","DECIPHER", "data.table") # "phangorn"

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)

dir <- '/Users/cigom/metagenomics/db/prep_model_DanioR'

setwd(dir)

file <- dir(pattern="derep.afa")

dna <- readDNAStringSet(file, format="fasta")

# Check intergenic noise
# freq <- letterFrequency(dna, letters="ACGTN-", OR=0)
# 
# library(superheat)
# superheat(freq,
#           scale = TRUE,
#           # add row dendrogram
#           row.dendrogram = TRUE)

#
# remove intergenic gapped sequences ----
# 

cleanSeqs <- function(DNAStringSet, pattern) {
  # it check and filter a set of intergenic insertions 
  
  pattern <- '[ATCG]N+[ATCG]|[ATCG]-+[ATCG]'
  noise <- DNAStringSet[grep(pattern = pattern, DNAStringSet)]
  
  rmNames <- names(noise)
  
  cleaned <- DNAStringSet[!names(DNAStringSet) %in% rmNames]
  
  return(cleaned)
}

# BrowseSeqs(noise)

cleaned_dna <- cleanSeqs(dna)

# BrowseSeqs(cleaned_dna)
# clean_freq <- letterFrequency(dna, letters="ACGTN-", OR=0)

outtag <- sub(".afa","",file)

writeXStringSet(cleaned_dna, 
                filepath = paste0(outtag, '_clean.afa'), 
                append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

# mandamos el alineamiento de secuencias 

quit(save = 'no' )


# Using decipher
# deshacemos el alinemiento

seqs <- RemoveGaps(dna)

plot(table(width(seqs)))

# name "Vertebrate Mitochondrial" ,
SGC1 <- getGeneticCode("SGC1")

gT <- lapply(order(width(seqs), decreasing=TRUE),
             function(x) {
               attr(x, "height") <- 0
               attr(x, "label") <- names(seqs)[x]
               attr(x, "members") <- 1L
               attr(x, "leaf") <- TRUE
               x
             })


attr(gT, "height") <- 0.5
attr(gT, "members") <- length(seqs)
class(gT) <- "dendrogram"
# use the guide tree as input for alignment

?AlignTranslation
DNA <- AlignTranslation(seqs,
                        guideTree=gT,
                        iterations=0,
                        refinements=0,
                        geneticCode = SGC1)


# DNA <- cleanSeqs(DNA)

#

BrowseSeqs(head(DNA))

writeXStringSet(DNA, 
                filepath = paste0(outtag, ".decipher.afa"), 
                append=FALSE,
                compress=FALSE, 
                compression_level=NA, format="fasta")


# Categorize inter-intra insertions


# it check and filter a set of intergenic insertions 
DNAStringSet <- dna

# 
# DNAStringSet2 <- dna[!names(dna) %in% names(DNAStringSet)]

categorizeMSA <- function(DNAStringSet, char = char) {
  
  require(DECIPHER)
  
  chars <- data.table(TerminalChar(DNAStringSet, char = char))
  chars$ncw <- width(RemoveGaps(DNAStringSet)) #Non-Chart-width
  
  chars[leadingChar + trailingChar >= 1, Categorize := 'Intra']
  
  chars[abs(difference - ncw) != 0  & leadingChar + trailingChar < 1, 
        Categorize := 'Inter']
  
  chars[abs(difference - ncw) == 0  & leadingChar + trailingChar == 0, 
        Categorize := 'none']
  
  chars[abs(difference - ncw) != 0  & leadingChar + trailingChar > 1, 
        Categorize := 'both']
  
  chars[, Expand := round(1 - (ncw / difference), digits = 3)]
  
  return(chars)
}

chars_mafft <- categorizeMSA(dna, char = '-')
table(chars_mafft$Categorize)

chars_decipher <- categorizeMSA(DNA, char = '-')
table(chars_decipher$Categorize)

library(ggplot2)
library(ggsci)

gg <- ggplot(chars_mafft, aes(difference, ncw, size = Expand, color = Categorize)) +
  geom_point() +
  theme_bw(base_size = 14) +
  scale_color_aaas() +
  labs(x = 'Sequence length (MSA)', 
       y = 'Sequence length (Barcodes)', 
       color = 'Group of\nInsertions',
       size = 'Contribution of\nInsertions',
       subtitle = paste0('Overview of'," '", char,"' ",'Isertions')) +
  theme(legend.position="bottom", legend.box = 'vertical') +
  facet_wrap(~ Categorize, ncol = 2)

tbl <- data.frame(table(chars_mafft$Categorize))
names(tbl) <- c('Category', 'Sequences')
library(patchwork)


box <- ggplot(chars_mafft, aes(y = ncw, x = )) +
  geom_boxplot() +
  theme_bw(base_size = 14) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.2) +
  labs(x = '', 
       y = 'Sequence length (Barcodes)')
  
  
  #stat_summary(fun.data=mean_sdl, mult=1, 
   #            geom="pointrange", color="red")

save <- gg + (box / gridExtra::tableGrob(tbl))

# gridExtra::tableGrob(mtcars[1:10, c('mpg', 'disp')])
ggsave(save, filename = paste0(dir, 'insertions.png'), width = 10, height = 7)

# E <- paste0('[A-Z]',char, '+[A-Z]') # Intergenic insertions
# A <-  paste0(char,'+[A-Z]') # Intragenic insertions


