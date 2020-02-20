
require(tidyverse)
require(ggplot2)
path <- 'metagenomics/train_rdp_classifier/fish_data/'

setwd(path)

file1 <- 'ROC_peces_disponibles_oksp_loso.txt'
file2 <- 'ROC_peces_disponibles_oksp_loto.txt'

# bootstrap	
# rank_FPR	
# rank_TPR	
# rank_F1score


# Calculate: 
# Accuracy (EXACTITUD) = VP + VN / VP + FP + FN + VN
# Precision = (TRP / TRP +  FPR)

rank <- c('Kingdom',	'Class',	'Order',	'Family',	'Genus',	'Species')
x <- read.table(file1, 
                sep="\t", skip = 2, header = T, 
                na.strings = "NaN")

names <- names(x)

paste0(rank, '_TPR') %in% names[grep('_TPR',names)]

get_TPR <- paste0(rank, '_TPR')
get_FPR <- paste0(rank, '_FPR')

x %>% select(get_TPR) %>% as_tibble() -> TPR
x %>% select(get_FPR) %>% as_tibble() -> FPR

TPR[is.na(TPR)] <- 0
FPR[is.na(FPR)] <- 0

# Precision = (TPR / TPR +  FPR)

P <- TPR / (TPR + FPR)
names(P) <- rank

png("boxplot.png", units="px", res=150)
boxplot(P, xlab, main = 'Precision\n(P = TPR / TPR + FPR)', las = 2)
dev.off()

plotLoot <- function(file) {
  
  options(stringsAsFactors = FALSE)
  
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  require(stringr)

  
  rank <- c('Kingdom',	'Class',	'Order',	'Family',	'Genus',	'Species')
  
  rScale <- c("Domain"="#edf8b1",  
              "Kingdom"="#7fcdbb", "Phylum"="#2c7fb8",  
              "Class"="#feb24c",  "Order"="#addd8e",  "Family"="#31a354",
              "Genus"="#bcbddc", "Species"="#756bb1")
  
  x <- read.table(file, 
                  sep="\t", skip = 2, header = T, 
                  na.strings = "NaN")
  
  names <- names(x)
  
  get_TPR <- names[grep('TPR',names)]
  get_FPR <- names[grep('FPR',names)]
  
  x %>% select(bootstrap, get_TPR) %>%
    melt(id.vars = 'bootstrap') %>%
    data.frame(ROC='True_Positive') -> TPR
  
  x %>% select(bootstrap, get_FPR) %>%
    melt(id.vars = 'bootstrap') %>%
    data.frame(ROC='False_Positive') -> FPR

  dtvs <- as_tibble(rbind(TPR, FPR))
  dtvs$Rank <- do.call(rbind, str_split(dtvs$variable, "_"))[,1]
  dtvs$Rank <- factor(dtvs$Rank, levels = rank)
  
  p <- ggplot(dtvs,
         aes(x=bootstrap, y=value, group=Rank, 
             fill=Rank, color=Rank)) + 
    geom_line(size=1, alpha=0.6, linetype=1) +
    scale_color_manual(values = rScale) +
    labs(title="", 
         caption = paste0("Results from file: ",file), 
         subtitle = "",
         x = 'Bootstrap',
         y = 'Rate') +
    theme_bw(base_size = 14) +
    facet_grid(~ ROC)
  
  filename <- str_split(file, "[.]")[[1]][1]
  
  ggsave(p, filename = paste0(filename, ".png"), dpi = 300,
         width = 7, height = 6)
  
  return(p)
}

plotLoot(file1)
plotLoot(file2)

# 



# library(caret)
# format my own data and run loot results
