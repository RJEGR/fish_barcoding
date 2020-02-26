devtools::install_github("TS404/AlignStat")

library("AlignStat")

data("reference_alignment")
data("comparison_alignment")

PAC <- compare_alignments(reference_alignment, 
                          comparison_alignment, CS=TRUE, SP=TRUE)
# Results visualisation
plot_similarity_heatmap    (PAC)
plot_dissimilarity_matrix  (PAC)
plot_similarity_summary    (PAC, CS=TRUE, cys=TRUE)
plot_dissimilarity_summary (PAC, stack=TRUE)
plot_SP_summary            (PAC, CS=TRUE)


# own
identical(names(dna), names(DNA))

mafft_align <- data.frame(t(as.matrix(dna)))
deciph_align <- data.frame(t(as.matrix(DNA)))

names(mafft_align) <- names(dna)
names(deciph_align) <- names(dna)

PAC <- compare_alignments(mafft_align, deciph_align, CS=TRUE, SP=TRUE) 

plot_similarity_heatmap    (PAC)
plot_dissimilarity_matrix  (PAC)
plot_similarity_summary    (PAC, CS=TRUE, cys=TRUE)
plot_dissimilarity_summary (PAC, stack=TRUE)
plot_SP_summary            (PAC, CS=TRUE)

