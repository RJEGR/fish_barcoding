path <- 'metagenomics/train_rdp_classifier/'

file1 <- 'metagenomics/train_rdp_classifier/ROC_Armatimonadetes_loso.txt'
file2 <- 'metagenomics/train_rdp_classifier/ROC_Armatimonadetes_loto.txt'

s <- read.table(file1, sep="\t", skip = 2, header = T, na.strings = "NaN")
t <- read.table(file2, sep="\t", skip = 2, header = T)

# bootstrap	
# rank_FPR	
# rank_TPR	
# rank_F1score
library(tidyverse)
library(reshape2)
library(ggplot2)

names <- names(s)

get_names <- names[grep('F1score',names)]

s %>% select(bootstrap, get_names) %>%
  melt(id.vars = 'bootstrap') %>%
  as_tibble() -> datavis
  

ggplot(datavis, 
       aes(x=bootstrap, y=value, group=variable, fill=variable, color=variable)) + 
  geom_line(size=1, alpha=0.6, linetype=1) +
  scale_color_brewer(palette = "Set1") +
  labs(title="",
       subtitle = '',
       x = 'Bootstrap',
       y = 'Accuracy') +
  theme_bw(base_size = 14)

# library(caret)
