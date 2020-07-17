# repasa este codigo. esta rarismo los numeros.
# hacer barras en base a taxa  y no a n secuencias!
# repasar el scaterplot de covertura
# elaborar sanity-check!

# Test groups of taxa 
# GROUP 1 ::
# Folmer Complete (FC) coverage, NCBI (N) and BOLD (B), 48, 705, 650
# GROUP 2 ::
# Folmer Partial (PC) coverage, NCBI (N), and BOLD (B), 48, 705, 500
# GROUP 3 ::
# Leray Complete 3p (LC) Coverage, NCBI (N) and BOLD B, 393, 705, 310
# GROUP 4 :: 
# Leray Partial 3p (LP) Coverage, NCBI (N) and BOLD (B), 393, 705, 260
# GROUP 5 :: 
# Leray 2 Complete 5p (JC) Coverage, NCBI (N) and BOLD (B), 48, 367, 310
# GROUP 6 :: 
# Leray 2 Partial 5p (JP) Coverage, NCBI (N) and BOLD (B), 48, 367, 310

# ----- analyze aligment
rm(list = ls())
options(stringsAsFactors = FALSE)

library(Biostrings)
library(tidyverse)
library(data.table)
library(DECIPHER)

# make function

summary_align <- function(DNAStringSet, start, end, range, factor) {
  
  rank <- c("Phylum","Class", "Order", "Family", "Genus", "Species")
  
  start <- start
  end <- end
  range <- range
  width <- end - start

  d <- subseq(DNAStringSet, start = start, end = end)

  # Determine the Number of Terminal Characters
  c_gaps <- data.table(DECIPHER::TerminalChar(d))
  
  seqs <- as.character(d)
  names(seqs) <- NULL

  d <- DECIPHER::RemoveGaps(d)
  tax <- data.frame(name = names(d))
  
  tax %>% separate(sep = ' ', 
                   col = 1, 
                   into = c('id', 'Taxonomy')) %>%
    separate(sep = ';',col = 2, into = rank) %>%
    as_tibble() -> tax
  
  cc <- data.table(start = start,
                   end = end,
                   width = width(d),
                   cov = width(d) / width,
                   int_cov = width(d) / c_gaps$difference,
                   names = tax$Species,
                   id = tax$id,
                   group = factor,
                   c_gaps,
                   seqs = seqs)
  
  cc <- cc[width >= range,]
  
  cc %>%
    as_tibble() %>%
    group_by(names) %>%
    group_size() -> size
  
  cc %>%
    as_tibble() %>%
    group_by(names) %>%
    group_keys() -> names
  
  summary <- data.table(size, names, group = factor)
  summary <- summary[order(-size), ]
  
  out <- list(summary = summary, res = cc)
  
  return(out)
  
}

# files ----
rank <- c("Phylum","Class", "Order", "Family", "Genus", "Species")
wr <- paste0(rank, "_wm")

rank_color <- c("Domain"="#edf8b1","Kingdom"="#7fcdbb", 
                "Phylum"="#2c7fb8","Class"="#feb24c",  
                "Order"="#addd8e",  "Family"="#31a354",
                "Genus"="#bcbddc", "Species"="#756bb1")

rf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher.afa'
bf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher_vs_metafishgom_bold.afa'
nf <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus.decipher_vs_metafishgom_ncbi.afa'

target <- 'metafishgom_DB_target.tsv'

p <- '/Users/cigom/metagenomics/db/fish_divergence/mafft/'
wd <- '/Users/cigom/metagenomics/db/fish_divergence/'
setwd(wd)
# Load files

ref <- readDNAStringSet(paste0(p, rf), format="fasta")
ncbi <- readDNAStringSet(paste0(p, nf), format="fasta") # 44947
bold <- readDNAStringSet(paste0(p, bf), format="fasta") # 38586

targetDB <- read.csv(target, sep = '\t')

dim(targetDB) # 7033

# morf <- c("morfologia_XIXIMI_01", "morfologia_XIXIMI_02", 
#          "morfologia_XIXIMI_03", "morfologia_XIXIMI_04",
#          "morfologia_XIXIMI_05", "morfologia_XIXIMI_06")

#morf_names <- paste0(substr(morf, 1,1),'x',
#                     substr(morf, nchar(morf),nchar(morf)))

pa <- c("bold_pa", "ncbi_pa", "BN_pa","morfo_pa")

targetDB %>%
  as_tibble() %>%
  mutate(Species_wm = str_replace(Species_wm, " ", "_")) %>%
  #rename(morf_names = morf) %>%
  select(AphiaID, wr, pa, GOM) -> gom

# names(gom)[names(gom) %in% morf] <- morf_names
# names(gom)[names(gom) %in% wr[3:6]] <- rank[3:6]

# Remove the reference from the aligment

bold <- bold[!names(bold) %in% names(ref)] # 38567
ncbi <- ncbi[!names(ncbi) %in% names(ref)] # 44928

# GROUP 1 ::
# Folmer Complete (FC) coverage, NCBI (N) and BOLD (B)
r <- 650
FCN <- summary_align(ncbi, 48, 705, r, 'FCN')
FCB <- summary_align(bold, 48, 705,  r, 'FCB')

# GROUP 2 ::
# Folmer Partial (PC) coverage, NCBI (N), and BOLD (B)
r <- 500
FPN <- summary_align(ncbi, 48, 705, r , 'FPN')
FPB <- summary_align(bold, 48, 705, r, 'FPB')

# GROUP 3 ::
# Leray 1 Complete (LC) Coverage, NCBI (N) and BOLD B
r <- 310
LCN <- summary_align(ncbi, 393, 705, r, 'LCN')
LCB <- summary_align(bold, 393, 705, r, 'LCB')

# GROUP 4 :: 
# Leray Partial (LP) Coverage, NCBI (N) and BOLD (B)
r <- 260
LPN <- summary_align(ncbi, 393, 705, r, 'LPN')
LPB <- summary_align(bold, 393, 705, r, 'LPB')

# GROUP 5 :: 
# Leray 2(L2) Coverage, NCBI (N) and BOLD (B)
r <- 315
L2N <- summary_align(ncbi, 48, 367, r, 'JCN')
L2B <- summary_align(bold, 48, 367, r, 'JCB')

# GROUP 6 :: 
# Leray 2 Partial (L2P) Coverage, NCBI (N) and BOLD (B)
r <- 260
L2PN <- summary_align(ncbi, 48, 367, r, 'JPN')
L2PB <- summary_align(bold, 48, 367, r, 'JPB')

# Test internal coverage
# res <- data.table(L2PN$res)
# res <- res %>% filter(int_cov < 1)
# 
# DNAset <- DNAStringSet(as.character(res$seqs))
# names(DNAset) <- res$id
# 
# DECIPHER::BrowseSeqs(DNAset)

summary <- rbind(
  FCN$summary,FCB$summary,
  FPN$summary, FPB$summary,
  LCN$summary, LCB$summary,
  LPN$summary, LPB$summary,
  L2N$summary, L2B$summary,
  L2PN$summary, L2PB$summary)

# Test summary 
n <- function(x) {length(unique(na.omit(x)))}

snf <- aggregate(summary[,'names'], by=list(summary$group), FUN = n)
sns <- aggregate(summary[,'size'], by=list(summary$group), FUN = sum)

# aggregate(summary$size, by = list(summary$group), function()

identical(snf[,1], sns[,1])

nsum <- data.table(group = snf[,1], sp = snf[,2], 
                   pct = snf[,2] / 7033, # este pct no es correcto!!
                   seqs = sns[,2])

nsum$set <- substr(nsum$group, 1,2)
nsum$marker <- substr(nsum$group, 1,1)
nsum[marker == 'L', marker:= 'Leray_3p']
nsum[marker == 'J', marker:= 'Leray_5p']
nsum[marker == 'F', marker:= 'Folmer']
#
nsum$db <- substr(nsum$group, nchar(nsum$group), nchar(nsum$group))
nsum[db == 'N', db:= 'GBank']
nsum[db == 'B', db:= 'BOLD']

write.table(nsum, file = 'summary_group.txt')


# group sequence sizes by ktone
summary[(size == 1 ), seq_group := "A"]
summary[(size > 1 & size <= 5), seq_group := "B"]
summary[(size > 5 & size <= 20), seq_group := "C"]
summary[(size > 20), seq_group := "D"]

# Cuantas especies tienen mas de una secuencia en cada grupo???:

ktf <- aggregate(summary[,'names'], by=list(summary$group, summary$seq_group), FUN = n)

kts <- aggregate(summary[,'size'], by=list(summary$group, summary$seq_group), FUN = sum)

identical(ktf[,1], kts[,1])

nktone <- data.table(group = ktf[,1], factor = ktf[,2], 
                             sp = ktf[,3] , 
                             seqs = kts[,3])

nktone$set <- substr(nktone$group, 1,2)
nktone$db <- substr(nktone$group, nchar(nktone$group), nchar(nktone$group))

nktone

nktone$marker <- substr(nktone$group, 1,1)
nktone[marker == 'F', marker:='Folmer']
nktone[marker == 'L', marker:='Leray_3p']
nktone[marker == 'J', marker:='Leray_5p']
nktone[marker == 'T', marker:='Target']

nktone[db == 'B', db:='BOLD']
nktone[db == 'N', db:='GBank']

summary[(size == 1 ), seq_group := "A"]
summary[(size > 1 & size <= 5), seq_group := "B"]
summary[(size > 5 & size <= 20), seq_group := "C"]
summary[(size > 20), seq_group := "D"]

labels <- c(A = "Singleton", B = "1 - 4", 
            C = "5 - 20", D = "> 20 seqs")

ggplot(nktone, aes(x = sp, y = seqs, 
                    color = set, shape = marker)) + 
  geom_point(size = 3, alpha = 0.7) + 
  #scale_color_brewer(type='div', palette='Set1') +
  scale_color_jco() +
  theme_bw() +
  theme(legend.position = "top") + 
  facet_grid(db ~ factor, scales = 'free_y', space="free",
             labeller = labeller(factor = labels)) + 
  labs(color = '', x = 'Number of species', y = 'Number of sequences')

ggsave("figure_2S.png", path = wd)

# Spread and save!

# spread <- select(nktone, -marker, -set, -factor,-sp)

# fish_encounters %>%
#   pivot_wider(names_from = station, values_from = seen)

#nktone %>% 
#  select(group, seqs) %>%
#  pivot_wider(names_from = group, values_from = seqs)

# Prepare data for bar-plot
summary <- select(summary, -seq_group)
# 0) check target groups
gom %>%
  filter(GOM == 1) %>%
  select(-Phylum_wm, -Class_wm) %>%
  melt(measure.vars = wr[3:6],
       variable.name = 'rank',
       value.name = 'name') %>% 
  data.table() -> target_only

target_only <- data.table(group =  'Target', target_only)

#

gom %>%
  select(-Phylum_wm, -Class_wm) %>%
  melt(measure.vars = wr[3:6],
       variable.name = 'rank',
       value.name = 'name') %>% 
  data.table() -> target_either

target_either <- data.table(group =  'Target', target_either)


# 1) check only reported species in gom

a <- filter(gom, GOM == 1)
  
right_join(summary, a, by = c('names' = 'Species_wm')) %>% 
  rename(Species_wm = names) %>% 
  select(-Phylum_wm, -Class_wm, -size) %>%
  melt(measure.vars = wr[3:6],
       variable.name = 'rank',
       value.name = 'name') %>% 
  data.table() -> only

nrow(only)

names(only) == names(target_only)
only <- rbind(only, target_only)

# 2) check either, gom and related species from genus reported into gom
b <- gom
left_join(summary, b, by = c('names' = 'Species_wm')) %>% 
  rename(Species_wm = names) %>%
  select(-Phylum_wm, -Class_wm, -size) %>%
  melt(measure.vars = wr[3:6],
       variable.name = 'rank',
       value.name = 'name') %>%
  data.table() -> either

names(either) == names(target_either)
either <- rbind(either, target_either)

# 3) check number of bold-ncbi related species
select(gom, -c(wr))

left_join(summary, b, by = c('names' = 'Species_wm')) %>% 
  rename(Species_wm = names) %>%
  select(-Phylum_wm, -Class_wm) %>% 
  filter(BN_pa == 1) %>%
  melt(measure.vars = wr[3:6],
       variable.name = 'rank',
       value.name = 'name') %>%
  data.table() -> presence_either

# 3.1 Presence only in gom

left_join(summary, b, by = c('names' = 'Species_wm')) %>% 
  rename(Species_wm = names) %>%
  select(-Phylum_wm, -Class_wm) %>% 
  filter(BN_pa == 1 & GOM == 1) %>%
  melt(measure.vars = wr[3:6],
       variable.name = 'rank',
       value.name = 'name') %>%
  data.table() -> presence_only

# Test groups of gom and either 
ncount_groups <- function(x) {
  n <- function(x) {length(unique(na.omit(x)))}

  ng <- aggregate(x[,'name'], by=list(x$rank, x$group), FUN = n)
  ng <- data.table(rank = ng[,1], 
                   group = ng[,2], ntaxa = ng[,3])
  
  ng$db <- substr(ng$group, nchar(ng$group),nchar(ng$group) )
  ng[db == 'N', db:='GBank']
  ng[db == 'B', db:='BOLD']
  
  ng$marker <- substr(ng$group, 1,1)
  ng[marker == 'F', marker:='Folmer']
  ng[marker == 'L', marker:='Leray_3p']
  ng[marker == 'J', marker:='Leray_5p']
  ng[marker == 'T', marker:='Target']
  
  ng$set <- substr(ng$group, 1,2)
  
  return(ng)
}

ng_a<- data.table(ncount_groups(only), wrap = 'Only')
ng_b <- data.table(ncount_groups(either), wrap = 'Either')

# special case for bold-ncbi set
ng_c <- data.table(ncount_groups(presence_either), wrap = 'Either')
ng_c <- aggregate(ng_c[,'ntaxa'], by=list(ng_c$rank, ng_c$set), FUN = sum)

marker <- data.table(marker = substr(ng_c[,2], 1,1))
marker[marker == 'F', marker:='Folmer']
marker[marker == 'L', marker:='Leray_3p']
marker[marker == 'J', marker:='Leray_5p']
marker[marker == 'T', marker:='Target']


ng_c <- data.table(rank = ng_c[,1], group = ng_c[,2],
                   ntaxa = ng_c[,3], db = 'bold-gbank', 
                   marker = marker$marker, set = ng_c[,2],
                   wrap = 'Either')

# 4)
ng_d <- data.table(ncount_groups(presence_only), wrap = 'Only')
ng_d <- aggregate(ng_d[,'ntaxa'], by=list(ng_d$rank, ng_d$set), FUN = sum)

marker <- data.table(marker = substr(ng_d[,2], 1,1))
marker[marker == 'F', marker:='Folmer']
marker[marker == 'L', marker:='Leray_3p']
marker[marker == 'J', marker:='Leray_5p']
marker[marker == 'T', marker:='Target']


ng_d <- data.table(rank = ng_d[,1], group = ng_d[,2],
                   ntaxa = ng_d[,3], db = 'bold-gbank', 
                   marker = marker$marker, set = ng_d[,2],
                   wrap = 'Only')


ng <- rbind(ng_a, ng_b, ng_c, ng_d)

write.table(ng, file = 'n_taxa_per_rank.csv', sep = '\t', row.names = FALSE)

# sanity check

ng$db <- factor(ng$db, levels = c('bold-gbank', 'BOLD', 'GBank'))

library(ggsci)

p1 <- ggplot(filter(ng, group != 'Target'), aes(x=rank, y=ntaxa, color = set, group = group, shape = marker)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_path() + 
  #facet_grid(db ~ wrap, scales = 'free_y', space = 'free_y', drop = FALSE) +
  facet_wrap(db~wrap, scales = 'free_y', ncol = 2)  +
  scale_color_jco() + 
  labs(title=" ",
       subtitle = '',
       color = NULL,
       x= 'Taxonomic level',
       y= "Number of taxa") + theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 <- p1 + scale_x_discrete(labels = c("Order_wm" = "Order",  
                                 "Family_wm"= "Family", 
                                 "Genus_wm" = "Genus",  
                                 "Species_wm" = "Species"))

ggsave("figure_1S.png", path = wd, width = 5.67,height =  7)

# Prepare barplot based on percent of target (either and only) groups
data.table(ncount_groups(target_only),  wrap = 'target')
data.table(ncount_groups(target_either), wrap = 'target')


ggplot(data = filter(ng, group == 'Target'), aes(x=rank, y=ntaxa, fill=set)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_jco()+
  theme_minimal(base_size = 14) + facet_grid(db~wrap, scales = 'free_y') +
  labs(title=" ",
       subtitle = '',
       fill = NULL,
       x= 'Taxonomic level',
       y= "Number of taxa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("Order_wm" = "Order",  
                              "Family_wm"= "Family", 
                              "Genus_wm" = "Genus",  
                              "Species_wm" = "Species"))

ggsave("figure_1S_bars.png", path = wd)

# results ----

res <- rbind(
  FCN$res,FCB$res,
  FPN$res, FPB$res,
  LCN$res, LCB$res,
  LPN$res, LPB$res,
  L2N$res, L2B$res,
  L2PN$res, L2PB$res)

write.table(FCN$res, file = 'res_example.csv', sep = '\t', row.names = FALSE)

# library(modeest) # carga la biblioteca modeest
# mfv(notas) # calcula la moda
# aggregate(res$cov, by = list(res$db, res$group), FUN = summary)

library(ggplot2)
# library(ggsci)
library(ggpubr)

# Visualize: Specify the comparisons you want

res$marker <- substr(res$group, 1,1)
res[marker == 'F', marker:='Folmer']
res[marker == 'L', marker:='Leray_3p']
res[marker == 'J', marker:='Leray_5p']
res[marker == 'T', marker:='Target']

res$db <- substr(res$group, nchar(res$group),nchar(res$group) )
res[db == 'N', db := "GBank"]
res[db == 'B', db := "BOLD"]

res$set <- substr(res$group, 1,2 )


g1 <- c("FCB", "FCN")
g2 <- c("FPB", "FPN")
g3 <- c("LCB", "LCN")
g4 <- c("LPB", "LPN")
g5 <- c("JCB", "JCN")
g6 <- c("JPB", "JPN")


my_comparisons <- list(g1, g2, g3, g4, g5, g6)


ggboxplot(res, x = "group", y = "cov",
          alpha = 0.5,
          group = 'group',
          color =  'set',
          palette = 'jco',
          ylab = 'Coverage (%)',
          xlab = 'Gene Group') + 
  stat_compare_means(comparisons = my_comparisons, 
                     label = "p.signif",
                     aes(color = set)) +
  stat_compare_means(label.y = 1.3)

ggsave("coverage_box.png", path = wd)


ggplot(res, aes(cov, fill = set, color = set)) + 
  #geom_histogram(alpha = 0.8, binwidth = 0.005, position = position_dodge2()) + 
  geom_freqpoly(alpha = 0.8, binwidth = 0.005, position = position_dodge2()) +
  scale_fill_jco() +
  scale_color_jco() +
  facet_grid(marker~db, scales = 'free_y') + theme_classic()

ggsave("coverage_freqpoly.png", path = wd)


# By number of seqs
summary$set <- substr(summary$group, 1,2)

ggboxplot(summary, x = "group", y = "size",
          alpha = 0.5,
          group = 'set',
          color =  'set',
          palette = 'jco',
          yscale = 'log2',
          notch = TRUE,
          ylab = '# Sequences (log2)',
          xlab = '') + 
  stat_compare_means(comparisons = my_comparisons, 
                     #label = "p.signif",
                     aes(color = set))

ggsave("sequence_box.png", path = wd)

x <- summary %>% 
  select(-set) %>% 
  spread(group, size, fill = 0)

rownames(x) <- x$names

x <- select(x, -names)

# subset to target db

db <- as_tibble(x, rownames = 'Species_wm')
targetDB %>%
  as_tibble() %>%
  mutate(Species_wm = str_replace(Species_wm, " ", "_")) %>%
  left_join(db, by = 'Species_wm') -> save

# save[is.na(save), ] <- 0

# save[save$GOM == 1, 'GOM'] <- 'T'
#save[GOM == 0, GOM := 'F']

write.table(save, file = 'metafishgom_DB_target_v2.tsv', sep = '\t', 
            row.names = FALSE)

library(reshape2)

left_join(gom, db, by = 'Species_wm') %>%
  melt(measure.vars = names(x),
       variable.name = 'group',
       value.name = 'seqs') %>% 
  melt(measure.vars = wr,
       variable.name = 'rank',
       value.name = 'name') %>% data.table() -> out

out$set <- substr(out$group, 1,2)
out$marker <- substr(out$group, 1,1)
out[marker == 'F', marker := "Folmer"]
out[marker == 'J', marker := "Leray_5p"]
out[marker == 'L', marker := "Leray_3p"]

# out[is.na(seqs), reads:= 0]

ggboxplot(out, x = "group", y = "seqs",
          alpha = 0.3,
          group = 'set',
          color =  'set',
          palette = 'jco',
          ylab = 'Number of Sequences',
          xlab = '') + 
  stat_compare_means(comparisons = my_comparisons, 
                     aes(color = set))


# metacoder! ----
library(viridis)
library(metacoder)
set.seed(20191201)

# Plasma
# magma

color_plasma <- c(viridis::plasma(3))
color_magma <- c(viridis::magma(3))

save[is.na(save)] <- 0

dim(input <- save %>% filter(GOM == 1))
samples <- unique(summary$group)

parse_tax_tbl <- function(tax_tbl, ranks, sam_data) {
  
  samples <- sam_data[[1]]
  
  tax_data <- as.data.frame(select(tax_tbl, ranks), 
                            stringsAsFactors = FALSE)
  
  tax_cols <- colnames(tax_data)
  tax_data <- cbind(data.frame(otu_id = rownames(tax_data), 
                               stringsAsFactors = FALSE), tax_data)
  
  otu_tbl <- as.data.frame(select(tax_tbl, samples), 
                           stringsAsFactors = FALSE)
  
  otu_tbl <- cbind(data.frame(otu_id = rownames(otu_tbl), 
                              stringsAsFactors = FALSE), otu_tbl)
  
  
  sam_data <- sam_data[sam_data[[1]] %in% samples]
  sam_data <- as.data.frame(sam_data, stringsAsFactors = FALSE)
  sam_data <- cbind(sample_id = rownames(sam_data), 
                    sam_data)
  
  datasets <- list()
  mappings <- c()
  
  datasets <- c(datasets, list(otu_table = otu_tbl))
  mappings <- c(mappings, c(`{{name}}` = "{{name}}"))
  
  datasets <- c(datasets, list(sample_data = sam_data))
  
  mappings <- c(mappings, NA)
  
  class_regex = "(.*)"
  class_key = "taxon_name"
  
  obj <- taxa::parse_tax_data(tax_data = tax_data, 
                              datasets = datasets, 
                              class_cols = tax_cols, 
                              mappings = mappings, 
                              named_by_rank = TRUE, 
                              #class_regex = class_regex,
                              class_key = class_key)
  
  return(obj)
}


obj <- parse_tax_tbl(input, wr, nsum)
# Getting per taxon information
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",
                                       cols = samples)

obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = samples,
                                      groups = sam_data$set)


# time demand!!

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs, 
                 node_label = taxon_names,
                 node_color = log2_median_ratio, 
                 node_color_range = diverging_palette(), 
                 node_color_trans = "linear", 
                 node_color_interval = c(-3, 3), 
                 edge_color_interval = c(-3, 3), 
                 node_size_axis_label = "Number of Species",
                 node_color_axis_label = "Log2 ratio median",
                 layout = "davidson-harel", 
                 initial_layout = "reingold-tilford",
                 output_file = "differential_heat_tree.pdf")

# puedes hacer porfa un metacoder solo con los Parciales? Que el color sea la presencia de secuencias y el tamano de nodo sea uniforme. A ver si queda un poquito mas clara la cobertura

samples <- nsum[[1]]
sam_partial <- samples[grep('P', samples)]
sam_data <- nsum[nsum[[1]] %in% sam_partial]
subset_in <- input[sample(nrow(input), 100), ]

obj <- parse_tax_tbl(input, wr, sam_data)

# Getting per taxon information
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",
                                       cols = sam_partial)

obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sam_partial,
                                      groups = sam_data$set)



# or!
obj %>%
  mutate_obs("tax_abund", 
             abundance = rowSums(obj$data$tax_abund[sam_data[["group"]]])) %>%
  heat_tree_matrix(data = "diff_table",
                 node_size = n_obs, 
                 node_label = taxon_names,
                 node_color = abundance, 
                 node_color_range = diverging_palette(), 
                 node_color_trans = "linear", 
                 node_size_axis_label = "Number of Species",
                 node_color_axis_label = "Number of Reads",
                 layout = "davidson-harel", 
                 initial_layout = "reingold-tilford",
                 output_file = "differential_heat_tree.pdf")


#
obj %>%
  filter_taxa(taxon_ranks == "Species_wm", supertaxa = TRUE) %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = root, 
            edge_color = leaf,
            initial_layout = "re", layout = "da",
            output_file = "plot_example.pdf")

# Primer set ----
# Folmer (650bp)

LCO1490 <- DNAString("GGTCAACAAATCATAAAGATATTGG")
HCO2198 <-  DNAString("TAAACTTCAGGGTGACCAAAAAATCA")

# Leray_2013 (313bp),
mlCOIintF <- DNAString("GGWACWGGWTGAACWGTWTAYCCYCC")
jgHCO2198 <- DNAString("TANACYTCNGGRTGNCCRAARAAYCA")

# match pattern!
table(res$set)

head(set <- filter(res, set == 'FC' & db  == 'BOLD'))
length(id <- set[['id']])
length(id <- unique(id))

start <- set[1,1]
end <- set[1,2]

dna <- bold

names <- unlist(lapply(strsplit(names(dna), 
                                      " ", fixed=TRUE), 
                             function(x) return(x[1])))
names(dna) <- names

dna <- subseq(dna, start, end)
dna <- dna[names(dna) %in% id]


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

# then run!!! the matchLeftAndRight!
coordAlign(mlCOIintF, dna[1])

