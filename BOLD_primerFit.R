# OBJETIVO CUANTAS ESPECIES ESTAN REPRESENTADAS EN EL FOLMER Y ESPECIES REPRESENTADAS EN EL FRAGMENTO LERAY, HOMOGENIZAR NOMENCLATURA A WORMS

rm(list = ls())

# task 
# evaluar tamanos de las secuencias y graficar en histograma, coloreando la ausencia/presencia total/parcial. Esto con las coordenadas de leray o folmer.

# alluvial / barplot de los casos a,b,c entre los primers folmer y leray (L1 y L2)


# Parameters ----
# (------------------------) profile, 658pb
# (------------------------) a, complete match
#             (------------) b, no match 5' , match in 3'
# (-----------)              c, match in 5' , no match in 3'
# ).                      .( d, no match 

# Define complete, partial and no match in left / right codon

# --- # ausencia total
# -NN # ausencia parcial
# --N # ausencia parcial
# NN- # ausencia parcial
# N-- # ausencia parcial
# NNN # presencia total

# for b and c, the limit of bp matchin will be 280, removiendo gaps despues de guardar las coordenas

# MAFFT INSTEAD OF HMMALIGN / HMMSCAN
# REMOVER SECUENCIAS CON GAPS
# - IDEALMENTE SE ESPERAN GAPS A LOS EXTREMOS (INTRA) Y NO INTERNAMENTE (INTER)
# - CALCULAMOS TRADUCCION DE A.A Y REMOVEMOS SECUENCIAS CON CODONES DE PARO INTERMEDIOS O AL INICIO
# - INNER_JOIN DE SECUENCIAS RESULTANTES E INFORMACION TAXONOMICA
# - ANALISIS DE DIVERGENCIA DE ESPECIES
# ETC.

# =======
# Library
# =======
.bioc_packages <- c("Biostrings", "ShortRead", "DECIPHER", "IRanges", "ggplot2") # "phangorn"

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)

#Inter-range features: ----
# input
#         (--------)
#              (----------------)
#                      (-----)
#                                     (---)
#                                     
# range   (-------------------------------)
# 
# reduce  (----------------------)    (---)
# 
# disjoin (---)(---)(-)(-----)(--)    (---)
# 
# GenomicRanges::setdiff(range(input),input)
#                                 (--)

# Intra-range features: ----
# input
#                         (----)
#                         .    .
# resize                  (--------)
# resize(fix="end")   (--------)
#                         .    .
# flank              (---).    .
# flank(start=F)          .    .(---)
#                         .    .
# promoters          (------)  .
#                         .    .
# narrow                  .(--).
#                         .    .
# shift (ignores strand!) .  (----)


# Bold model alignment ----

# with hhm
# c('/Users/cigom/metagenomics/db/bos_taurus_align/rep_otu_bos_taurus/coi_profile_vs_peces_bold_no_gaps.hmm.align')

f1 <- c('/Users/cigom/metagenomics/db/mafft/peces_bold_0.2_cntrds_100_COI_bos_taurus_clean_vs_peces_bold.afa')

f2 <- c('/Users/cigom/metagenomics/db/mafft/mafft_clean_vs_peces_bold.afa')

f3 <- c('/Users/cigom/metagenomics/db/bos_taurus_align/rep_otu_bos_taurus/coi_profile_vs_peces_bold_no_gaps.hmm.align')

# file <- dir(path = '/Users/cigom/metagenomics/db/bos_taurus_align/coi_refaln/', pattern = "align", full.names = TRUE)


seqs1 <- readDNAStringSet(f1, format="fasta")
seqs2 <- readDNAStringSet(f2, format="fasta")
seqs3 <- readDNAStringSet(f3, format="fasta")


library(data.table)

s1 <- data.table(seqSize(f1, 367, 731), f = '0.2_cntrds_100')
s2 <- data.table(seqSize(f2, 367, 731), f = 'mafft_clean')
s3 <- data.table(seqSize(f3, 367, 731), f = 'hmm')

sizes <- rbind(s1, s2, s3)
names(sizes) <- c('Width', 'Seqs', 'f')

sizes[, Width := as.numeric(Width)]
sizes[, Seqs := as.numeric(Seqs)]


sizes %>%
  ggplot( aes(x=Width, y=Seqs, group=f, fill=f)) +
  geom_area() + 
  scale_fill_viridis_d() + theme_bw(base_size = 14) +
  theme(legend.position="top") # facet_wrap(~f, scale="free_y") 


# pos-processing step 
# ref http://www2.decipher.codes/Documentation/ArtOfAlignmentInR.pdf
# with the goal of removing artifacts of the progressive alignment process. 
# This function will efficiently correct most obvious inaccuracies that could be found by-eye. 

seqSize <- function(fasta_path, start, end) {
  
  seqs <- readDNAStringSet(fasta_path, format="fasta")
  # seqs <- AdjustAlignment(seqs)
  coor <- subseq(seqs, start = start, end = end)
  coor <- DECIPHER::RemoveGaps(coor)
  
  ungapp <- DECIPHER::RemoveGaps(coor)
  size_tbl <- data.frame(table(width(ungapp)))
  
  return(size_tbl)
  
}

#Fragmento / primer Leray_2013 (313bp)
#* Inicio 393 o 394 / 367 (+26)
#* Termino 705 o 702 / 731 (+26)


# Primers Folmer_1994 ----
# LCO1490 long 25nt Inicio pos-23, termina pos-47
# HCO2198 (RevComp) long 26nt Inicio pos-706, termina pos-731

folmer <- subseq(seqs_adj, start = 23, end = 731)
folmer_ungapp <- DECIPHER::RemoveGaps(folmer) # si removemos gaps

barplot(table(width(folmer_ungapp)))


barplot(table(width(leray)))

BrowseSeqs(sample(folmer, 100), highlight=1)

# remove sequence less thann 280
widths <- width(leray)
sum(widths >= 280)
leray <- leray[widths >= 280]

BrowseSeqs(sample(leray, 100), highlight=0)

Lcodon <- subseq(leray, start = 1, end = 3)

lengths <- width(leray)
Rcodon <- subseq(leray, lengths - 2, lengths)

Ltb <- data.frame(table(Lcodon))
Rtb <- data.frame(table(Rcodon))

Ltb$Region <- 'Left'
Rtb$Region <- 'Rigth'

names(Rtb)[1] <- 'Oligo'
names(Ltb)[1] <-  'Oligo' 

tb <- rbind(Ltb, Rtb)

# --- # ausencia total
# -NN # ausencia parcial
# --N # ausencia parcial
# NN- # ausencia parcial
# N-- # ausencia parcial
# NNN # presencia total


ggplot(tb, aes(x = Oligo,y = Freq, fill = Region)) + geom_bar(stat = "identity", position = "stack", color = "black") + coord_flip() + theme_bw() + labs(y = 'Number of sequences', 'Oligo') + facet_grid(rows = vars(Region), space = 'free', scales = 'free') +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.1, hjust = -0.3)

ggplot(tb, aes(x=Oligo, y=Freq, size = Freq, color = Region)) +
  geom_point(alpha=0.5) + coord_flip() + theme_bw() + labs(y = 'Number of sequences', 'Oligo') +
  #facet_grid(rows = vars(Region), space = 'free', scales = 'free') + 
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), size = 3, vjust=0.1, hjust = 0.9) + ggtitle('Number of Oligos at start and end position in Leray fragment coordinates')

# Prepare bos taurus sequence using differnt coords. make a input of 10 seqs

seqs <- readDNAStringSet(bhmm_name, format="fasta")

bos <- readDNAStringSet("~/metagenomics/db/bos_taurus_align/coi_refaln/COI_bos_taurus.fasta")


# make different fragmetns of the bos based on the width and a normal distribution of 10 sizes
1:width(bos)

b1 <- subseq(bos, 1, 100)
b1 <- subseq(bos, 50, 100)


