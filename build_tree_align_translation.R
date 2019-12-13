# build_tree_align_translation.R 
# November 2019
# Ricardo Gomez | rgomez@cicese.mx

# =======
# Library
# =======

rm(list = ls())

.bioc_packages <- c("Biostrings","DECIPHER", "IRanges") # "phangorn"

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)

# refine aligment
# Aligning hundreds of thousands of unique sequences?

# Chained guide trees offer a viable
# alternative when aligning hundreds
# of thousands of unique sequences,
# as the default requires O(n2) time.
# Refer to the section of this vignette
# entitled “Building a Guide Tree”

# Load inputs ----
opath <- getwd()

args = commandArgs(trailingOnly=TRUE)

dna.file <- args[1]
seqs.file <- args[2]

# dna.file <- c("~/metagenomics/db/bos_taurus_align/coi_refaln/coi_btarurus_ictioconsenso_primers_aln.fasta")
seqs.file <- c("~/metagenomics/db/prep_model_ref/tmp_files_peces_bold/peces_bold_0.2_cntrds.tmp")
bos_taurus <- c("~/metagenomics/db/prep_model_ref/COI_bos_taurus.fasta")
# form a chained guide tree ----

dna <- readDNAStringSet(c(bos_taurus, seqs.file), format="fasta")
dna <- RemoveGaps(dna)
# dna <- sample(dna, 1000)
dna

gT <- lapply(order(width(dna), decreasing=TRUE),
             function(x) {
               attr(x, "height") <- 0
               attr(x, "label") <- names(dna)[x]
               attr(x, "members") <- 1L
               attr(x, "leaf") <- TRUE
               x
             })


attr(gT, "height") <- 0.5
attr(gT, "members") <- length(dna)
class(gT) <- "dendrogram"

# use the guide tree as input for alignment


DNA <- AlignTranslation(dna,
                        guideTree=gT,
                        iterations=0,
                        refinements=0)

#BrowseSeqs(DNA)

writeXStringSet(DNA, filepath = paste0(opath,"/", seqs.file, ".tree"), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

quit(save = 'no')

# 

folmer <- subseq(DNA, start = 23, end = 731)
folmer <- RemoveGaps(folmer)

BrowseSeqs(folmer)







