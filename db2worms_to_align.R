
# 1) Generar un formato fasta con la nueva nomenclatura worms, 
# 2) entonces hacemos un archivo filtrado de la lista de especies del gom

rm(list = ls())

options(stringsAsFactors = FALSE)

dir <- '/Users/cigom/metagenomics/db/fish_divergence'
setwd(dir)

qr <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

wr <- paste0(qr, "_wm")

# qr_ncbi <- c("phylum", "class", "order", "family", "genus", "Species")

bold2w <- read.csv('metafishgom_B2W.tsv', sep = '\t')
ncbi2w <- read.csv('metafishgom_subN2W.tsv', sep = '\t')

library(dplyr)
library(stringr)
library(tidyr)

bold2w %>%
  as_tibble() %>%
  select(processid, wr) %>% 
  mutate(Species_wm = str_replace(Species_wm, " ", "_")) %>%
  unite('Taxonomy', wr, sep = ';') -> b2w

# 2

ncbi2w %>%
  as_tibble() %>%
  select(accnum, wr) %>% 
  mutate(Species = str_replace(Species_wm, " ", "_")) %>%
  unite('Taxonomy', wr, sep = ';') -> n2w

# Parse sequence to taxonomy ----

ncbiF <- 'ncbi_complete.fasta'
boldF <- 'peces_bold_sorted.tmp'

boldL <- 'metafishgom_bold_processid.list'
ncbiL <- 'metafishgom_ncbi_accnum.list'

nrow(boldL <- read.table(boldL)) # 38, 575
nrow(ncbiL <- read.table(ncbiL)) # 44, 928

filter_fa <- function(fasta, list) {
  
  library('Biostrings')
  
  db <- readDNAStringSet(fasta, format="fasta")
  db <- db[names(db) %in% list[,1]]
  return(db)
}

ncbi <- filter_fa(ncbiF, ncbiL)
bold <- filter_fa(boldF, boldL)

# sanity check

n2w <- n2w[match(names(ncbi), n2w$accnum),]
b2w <- b2w[match(names(bold), b2w$processid),]

identical(names(ncbi), n2w$accnum)
identical(names(bold), b2w$processid)

# Save files

nn <- unite(n2w, 'names', accnum:Taxonomy, sep = ' ')
bn <- unite(b2w, 'names', processid:Taxonomy, sep = ' ')
  
names(ncbi) <- nn$names
names(bold) <- bn$names

writeXStringSet(bold, filepath = paste0(dir, '/metafishgom_bold.fa'), 
                append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


writeXStringSet(ncbi, filepath = paste0(dir, '/metafishgom_ncbi.fa'), 
                append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

# Create tiles by group ----
# if group by Order-groups?
ncbi2w %>%
  as_tibble() %>%
  select(accnum, wr) %>%
  group_by(Order_wm) %>%
  select(Order_wm) %>% 
  table() %>%
  data.table() -> ncbi2w_group

ncbi2w_group[order(N)]

bold2w %>%
  as_tibble() %>%
  select(processid, wr) %>%
  group_by(Order_wm) %>%
  select(Order_wm) %>% 
  table() %>%
  data.table() -> bold2w_group

bold2w_group[order(N)]

# n_groups(ncbi2w_group) # 

library(data.table)

x <- data.table(db = 'GB', table(ncbi2w_group[,1]))
y <- data.table(db = 'BOLD', table(bold2w_group[,1]))

GB <- data.table(db = 'GB', 
           group_keys(ncbi2w_group),
           size = group_size(ncbi2w_group))

BOLD <- data.table(db = 'BOLD', 
                 group_keys(bold2w_group),
                 size = group_size(bold2w_group))



db <- rbind(GB,BOLD, use.names=FALSE)
db <- db %>% spread(db, size, fill = NA)


#pct <- function(x) {x / sum(x) * 100}
#db[, BOLD := (pct(BOLD))]
#db[, GB := (pct(GB))]

# sanity check

colSums(db[,-1])

library(ggalluvial)
# 
# ggplot(data = db,
#        aes(axis1 = BOLD, axis2 = order ,
#            axis3 = GB)) +
#   geom_alluvium() +
#   scale_x_discrete(limits = c("BOLD", "Group", "Gene Bank")) +
#   geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
#   theme_minimal(base_size=16) +
#   labs(title = "", x = "" , y = "n_groups")


quit(save = 'no')