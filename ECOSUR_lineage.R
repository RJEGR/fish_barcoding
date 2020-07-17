#
# En estas secuencias se tiene que revisar si hay taxones nuevos, 
# es decir, que no esten en (la nueva version de) DB_target

# vamos:
# 1) contar niveles repetidos
# 2) comprar linajes en comun por metodo:
# 2.1) inner_join
# 2.2) Hamming Compute distance metrics between strings of taxa nomenclature

options(stringsAsFactors = FALSE)

path <- '/Users/cigom/metagenomics/db/BOLD_ECOSUR'

setwd(path)

# for i in $(ls *fasta); do grep '^>' $i| sed 's/>//g'  > ${i%.fasta}.lineage; done

# bold_ecosur_ISABZ.lineage
# bold_ecosur_MFLS.lineage
# bold_ecosur_MFLS.lineage

rtb <-function(file) {
  
  group <- str_split(file, "[.]")[[1]][1]
  group <- str_split(group, "[_]")[[1]][3]
  
  x<- read.csv(file, sep = "|", header = F, 
               fill=T, stringsAsFactors=FALSE)
  x<- data.frame(x, group = group)
}

x <- lapply(dir(pattern = 'lineage'), rtb)
x <- do.call(rbind, x)

names <- c('Sampleid', 'Bilateral','Processid', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Lifestage', 'Vouchertype', 'Collectors', 'Collectiondate', 'Country', 'Province', 'Region', 'Site', 'Lat', 'Long', 'Depth', 'Gps_source', 'Habitat', 'Sampling_protocol', 'COI-5P', 'BIN', 'URI', 'group')

rScale <- c("Domain"="#edf8b1",  
            "Kingdom"="#7fcdbb", "Phylum"="#2c7fb8",  
            "Class"="#feb24c",  "Order"="#addd8e",  
            "Family"="#31a354",
            "Genus"="#bcbddc", "Species"="#756bb1")

names(x) <- names

table(x$group)

ranks <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

require(reshape2)

x %>% select(ranks, group) %>% as_tibble()-> tax 

tax %>% 
  melt(id.vars = 'group', 
       value.name = 'lineage', 
       variable.name = 'rank') %>%
  as_tibble() -> shaped

# 2)
n <- function(x) {table(unique(x))}

agg_shaped <- as_tibble(aggregate(shaped$rank, 
          by = list(shaped$group,shaped$rank), 
          FUN = length ))

names(agg_shaped) <- c('group', 'rank', 'n') #'lineage',
  

ggplot(agg_shaped, 
      aes(x=rank, y=n, group=group, fill=group, color=group)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") + 
  theme_bw(base_size = 14)

# 

library(ggalluvial)

alluv <- to_lodes_form(tax, axis = 1:6)

# data(majors)
# head(majors)
# ggplot(majors, aes(x = semester, stratum = curriculum, alluvium = student, fill = curriculum))

ggplot(alluv,
       aes(x = group, stratum = x, 
           alluvium = alluvium, y = 1)) +
  geom_alluvium(fill = "darkgrey", na.rm = TRUE) +
  geom_stratum(aes(fill = x), color = NA, na.rm = TRUE) +
  theme_bw() +
  geom_flow(stat = "alluvium", color = "black") +
  scale_y_reverse() +
  scale_fill_manual(values=rScale)

# verificar como mejorar
# http://corybrunson.github.io/ggalluvial/reference/stat_alluvium.html
