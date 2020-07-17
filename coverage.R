
library(Biostrings)
library(tidyverse)

# options(stringAsFactor = FALSE)

fasta <- dir(pattern='fasta$', path='~/metagenomics/metafishgom_mocks', full.names=TRUE)

fl.bam <- dir(pattern=".bam$", path='~/metagenomics/metafishgom_mocks/coverage/bam', full.names=TRUE)

wd <- '~/metagenomics/metafishgom_mocks/coverage'

#ref <- readDNAStringSet(fastaFile)
#id <- sapply(strsplit(names(ref), " "), `[`, 1)

#widths <- width(ref)
#GC <- rowSums(letterFrequency(ref, letters="GC", OR=0)) / widths

# wh <- GRanges(seqnames = Rle(id),
#               ranges = IRanges(1, end = widths, names = id),
#               Rle(strand(c("+"))),
#               GC = GC)

library(Rsamtools)

ConBam <- fl.bam[grep('Con', fl.bam)]
DivBam <- fl.bam[grep('Div', fl.bam)]
RelBam <- fl.bam[grep('Rel', fl.bam)]

# bamList <- lapply(DivBam, function(x) scanBam(x))

bam2tbl <- function(bamfile) {
  
  bam <- scanBam(bamfile)
  
  .unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
      # structure(unlist(x), class = "character")
    } else {
      do.call(c, x)
    }
  }
  
  bam_field <- names(bam[[1]])
  
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field

  
  return(bam_df)
  
}

# ConBam rep ----
# con_df1 <- bam2tbl(ConBam[1])
# con_df2 <- bam2tbl(ConBam[2])
# con_df3 <- bam2tbl(ConBam[3])

Con_df <- rbind(data.frame(bam2tbl(ConBam[1]), Rep = 'rep1'),
                data.frame(bam2tbl(ConBam[2]), Rep = 'rep2'),
                data.frame(bam2tbl(ConBam[3]), Rep = 'rep3'))
# DivBam rep ----
# Div_df1 <- bam2tbl(DivBam[1])
# Div_df2 <- bam2tbl(DivBam[2])
# Div_df3 <- bam2tbl(DivBam[3])

Div_df <- rbind(data.frame(bam2tbl(DivBam[1]), Rep = 'rep1'),
                data.frame(bam2tbl(DivBam[2]), Rep = 'rep2'),
                data.frame(bam2tbl(DivBam[3]), Rep = 'rep3'))
# RelBam rep ----
# Rel_df1 <- bam2tbl(RelBam[1])
# Rel_df2 <- bam2tbl(RelBam[2])
# Rel_df3 <- bam2tbl(RelBam[3])

Rel_df <- rbind(data.frame(bam2tbl(RelBam[1]), Rep = 'rep1'),
                data.frame(bam2tbl(RelBam[2]), Rep = 'rep2'),
                data.frame(bam2tbl(RelBam[3]), Rep = 'rep3'))


# check refseq cov ----
data.frame(table(Con_df$rname))
data.frame(table(Div_df$rname))
data.frame(table(Rel_df$rname))
# flags
data.frame(table(Con_df$flag), Mock = 'Con')
data.frame(table(Div_df$flag), Mock = 'Div')
data.frame(table(Rel_df$flag), Mock = 'Rel')

strand_mapped <- function(rname, bam_df) {
  
  # decipher aligment information by check_* function
  #function for checking negative strand
  check_neg <- function(x){
    if (intToBits(x)[5] == 1){
      return(T)
    } else {
      return(F)
    }
  }
  
  #function for checking positive strand
  check_pos <- function(x){
    if (intToBits(x)[3] == 1){
      return(F)
    } else if (intToBits(x)[5] != 1){
      return(T)
    } else {
      return(F)
    }
  }
  
  #store the mapped positions on the plus and minus strands ====
  
  # bam_df %>%
  #   as_tibble() %>%
  #   mutate(flag = apply(as.data.frame(flag), 1, check_pos)) %>%
  #   filter(rname == rname) %>%
  #   mutate(rname = fct_drop(rname)) %>%
  #   select(pos, Rep, strand) -> rname_pos
  # 
  # bam_df %>%
  #   as_tibble() %>%
  #   mutate(flag = apply(as.data.frame(flag), 1, check_neg)) %>%
  #   filter(rname == rname & flag == TRUE) %>%
  #   mutate(rname = fct_drop(rname)) %>%
  #   select(pos, Rep, strand) -> rname_neg
  
  # previos version -----
  # rname_neg <- bam_df[bam_df$rname == rname &
  #                         apply(as.data.frame(bam_df$flag), 1, check_neg),
  #                       c('pos','Rep', 'strand')]
  # 
  # rname_pos <- bam_df[bam_df$rname == rname &
  #                       apply(as.data.frame(bam_df$flag), 1, check_pos),
  #                     c('pos','Rep', 'strand') ]
  # continue ----
  
  bam_df %>%
    filter(rname == rname) %>%
    mutate(rname = fct_drop(rname)) -> bam_df
  
  neg <- bam_df[apply(as.data.frame(bam_df$flag), 1, check_neg), 'pos']
  pos <- bam_df[apply(as.data.frame(bam_df$flag), 1, check_pos), 'pos'] 

  neg_density <- density(neg)
  pos_density <- density(pos)
  
  #display the negative strand with negative values
  neg_density$y <- neg_density$y * -1
  
  # return (list(neg_density, pos_density))
  
  # 
  # rep_neg <- as.character(rname_neg$Rep)
  # rep_pos <- as.character(rname_pos$Rep)
  # 
  # strand_neg <- as.character(rname_neg$strand)
  # strand_pos <- as.character(rname_pos$strand)
  # 
  
  # Rep <- as.character(unique(bam_df$Rep))
  
  yneg <- neg_density$y
  xneg <- neg_density$x
  
  ypos <- pos_density$y
  xpos <- pos_density$x
  
  # tbl <- rbind(data.frame(strand = "+", pos = rname_pos$pos, 
  #                         rname = rname, Rep = rep_pos),
  #              data.frame(strand = "-", pos = rname_neg$pos, 
  #                         rname = rname, Rep = rep_neg)) 
  
  tbl <- rbind(data.frame(y = ypos, x = xpos),
               data.frame(y = yneg, x = xneg))
  
  tbl <- data.frame(tbl, rname)

  return(tbl)
}


bam_df <- Con_df

rnames <- table(bam_df$rname) # as factor
rnames <- names(rnames[rnames > 0])
rname <-  rnames[12] # 'R029'


strand_mapped(rname = rname, bam_df) %>% as_tibble() %>%
  ggplot(aes(x,y, fill = rname)) +
  geom_area()

# collapse and plot ----

strand_tbl <- function(bam_df) {
  

  rnames <- table(bam_df$rname) # as factor
  rnames <- names(rnames[rnames > 0]) # if is factor , clean 
  
  bam_df %>%
    as_tibble() %>%
    filter(rname %in% rnames) %>%
    mutate(rname = fct_drop(rname)) -> bam_df
  
  allframes <- lapply(rnames, function(x) {strand_mapped(rname = x, bam_df)})
  
  tbl <- do.call(rbind, allframes) %>% as_tibble()
  
  return(tbl)
  
}

tbl_con <- strand_tbl(Con_df)

tbl_div <- strand_tbl(Div_df)

tbl_rel <- strand_tbl(Rel_df)



library(ggridges)

p1 <- tbl_con %>%
  ggplot(aes(x = pos, y = rname, fill = strand)) + # fill  = 
  geom_density_ridges(alpha = 0.5) +
  theme_ridges(font_size = 7) +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  labs(y = 'RefSeqs', x = ' RefSeqs Position')
  #facet_grid(~Rep)

p2 <- tbl_div %>%
  ggplot(aes(x = pos, y = rname, fill = strand)) + # fill  = 
  geom_density_ridges(alpha = 0.5) +
  theme_ridges(font_size = 7) +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  labs(y = 'RefSeqs', x = ' RefSeqs Position')

p3 <- tbl_rel %>%
  ggplot(aes(x = pos, y = rname, fill = strand)) + # fill  = 
  geom_density_ridges(alpha = 0.5) +
  theme_ridges(font_size = 7) +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  labs(y = 'RefSeqs', x = ' RefSeqs Position')

ggsave(p1, filename = 'coverage_con.png', path = wd, width = 10, height = 5.89)
ggsave(p2, filename = 'coverage_div.png', path = wd, width = 10, height = 5.89)
ggsave(p3, filename = 'coverage_rel.png', path = wd, width = 10, height = 5.89)



# tbl_con %>%
#   ggplot(aes(x = pos, y = Rep, fill = factor(stat(quantile)))) +
#   stat_density_ridges(
#     geom = "density_ridges_gradient",
#     calc_ecdf = TRUE,
#     quantiles = c(0.025, 0.975)
#     ) +
#   scale_fill_manual(
#     name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
#     labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]"))

# con lo siguiente te das cuenta que la replicas tienen stats similares

#aggregate(con_df$pos, by = list(con_df$Rep, con_df$rname), summary)

# using other approach ----
#https://stackoverflow.com/questions/53699250/merging-multiple-lists-into-a-data-frame-for-ggpot2
extracting_pos_neg_reads <- function(bam_fn) {
  
  #read in entire BAM file
  bam <- scanBam(bam_fn)

  #function for collapsing the list of lists into a single list
  #as per the Rsamtools vignette
  .unlist <- function (x) {
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)) {
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  #store names of BAM fields
  bam_field <- names(bam[[1]])
  
  #go through each BAM field and unlist
  list <- lapply(bam_field, function(y)
    .unlist(lapply(bam, "[[", y)))
  
  #store as data frame
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  
  #function for checking negative strand
  check_neg <- function(x) {
    if (intToBits(x)[5] == 1) {
      return(T)
    } else {
      return(F)
    }
  }
  
  #function for checking positive strand
  check_pos <- function(x) {
    if (intToBits(x)[3] == 1) {
      return(F)
    } else if (intToBits(x)[5] != 1) {
      return(T)
    } else {
      return(F)
    }
  }

  
  #store the mapped positions on the plus and minus strands
  neg <- bam_df[apply(as.data.frame(bam_df$flag), 1, check_neg),'pos']
  pos <- bam_df[apply(as.data.frame(bam_df$flag), 1, check_pos),'pos']
  
  #calculate the densities
  neg_density <- density(neg)
  pos_density <- density(pos)
  
  #display the negative strand with negative values
  neg_density$y <- neg_density$y * -1
  
  return (list(neg_density, pos_density))
  
}

all_densities <- data.frame()

groups <- c(paste0('Con', 1:3),
            paste0('Div', 1:3),
            paste0('Rel', 1:3))

k <- 1

for (i in fl.bam){ 
  print(i)
  #density <- extracting_pos_neg_reads(i)
  densities <- cbind(rbind(data.frame(density[[1]][1:2]), data.frame(density[[2]][1:2])), 
                     id = rep(c("neg", "pos"), each = length(density[[1]]$x)))
  densities$group <- groups[k]
  k <- k + 1
  all_densities <- rbind(all_densities, densities)
}

all_densities %>%
  filter(group %in% paste0('Con', 1:3)) %>%
  as_tibble() %>%
  ggplot(aes(x,y, fill = id)) +
  geom_area(alpha=0.4) +
  theme_bw() -> p1

# 2
all_densities %>%
  filter(group %in% paste0('Div', 1:3)) %>%
  as_tibble() %>%
  ggplot(aes(x,y, fill = id)) +
  geom_area(alpha=0.4) +
  # geom_line(aes(color=id), size=0.4) +
  # geom_point(size=0.4, aes(color = id)) +
  theme_bw() -> p2

all_densities %>%
  #filter(group %in% paste0('Rel', 1:3)) %>%
  as_tibble() %>%
  ggplot(aes(x,y, fill = id)) +
  geom_area(alpha=0.4, position = position_dodge()) +
  # geom_line(aes(color=id), size=0.4) +
  # geom_point(size=0.4, aes(color = id)) +
  theme_bw() -> p

ggsave(p, filename = 'coverage_mocks.png', path = wd, width = 10, height = 5.89)


