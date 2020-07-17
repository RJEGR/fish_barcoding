# devtools::install_github("TS404/AlignStat")
# 
# library("AlignStat")
# 
# data("reference_alignment")
# data("comparison_alignment")
# 
# PAC <- compare_alignments(reference_alignment, 
#                           comparison_alignment, CS=TRUE, SP=TRUE)
# # Results visualisation
# plot_similarity_heatmap    (PAC)
# plot_dissimilarity_matrix  (PAC)
# plot_similarity_summary    (PAC, CS=TRUE, cys=TRUE)
# plot_dissimilarity_summary (PAC, stack=TRUE)
# plot_SP_summary            (PAC, CS=TRUE)
# 

# own ----
library("AlignStat")

paste0(outtag, ".decipher.afa")

identical(names(dna), names(DNA))

mafft_align <- data.frame(t(as.matrix(dna)),  stringsAsFactors=TRUE)
deciph_align <- data.frame(t(as.matrix(DNA)), stringsAsFactors=TRUE)

#

names(mafft_align) <- 1:length(names(mafft_align))
names(deciph_align) <- 1:length(names(deciph_align))

PAC <- compare_alignments(head(mafft_align[1:1000], 10), 
                          head(deciph_align[1:1000], 10), CS=TRUE, SP=TRUE)

plot_similarity_heatmap    (PAC)
plot_dissimilarity_matrix  (PAC)
plot_similarity_summary    (PAC, CS=TRUE, cys=TRUE)
plot_dissimilarity_summary (PAC, stack=TRUE)
plot_SP_summary            (PAC, CS=TRUE)

