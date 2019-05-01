# Script for dealing with newly annotated Euk OTUs from 18S dataset
# Cariaco cruises 212 and 216
# New analysis started 2018-07-02

# Requires tidyverse, ggthemes (ggplot is in tidyverse) and vegan. First time install:
# install.packages("tidyverse") 
# install.packages("vegan") 
# install.packages("ggthemes")
# install.packages("devtools")
# install.packages("ggpubr")
# install.packages("reshape2")

# rm(list = ls())

# Load packages 
library(tidyverse)
library(vegan)
library(ggthemes)
library(gridExtra)
library(grid)
#devtools::install_github("baptiste/egg")
library(egg)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(cowplot)
library("igraph") 
library("network") 
library("sna") 
library("ndtv")
library("RColorBrewer") 
# display.brewer.all()
library(extrafont)
## Import system fonts - may take a while
# font_import() 
# fonts() 

# Read in new PR2 annotations
PR2_annotation <- read_delim("PR2_annotation/Cariaco_AE_qiime_updated_long.fa_rep_set_tax_assignments.txt",delim = "\t", col_names = T)

# Read in abundance table with old SILVA annotations
abun_table <- read_delim("iTags_Analyses_SILVA_old/Cariaco_AE_updated_long_sums.txt",delim = "\t", col_names = T)

# Delete SILVA annotations
abun_table <- abun_table[,1:49]

# Match PR2 annotations to OTU table
otu.table <- PR2_annotation %>% 
  left_join(abun_table, by = "#OTU ID") 

# separate taxonomy into multiple columns
otu.table <- separate(otu.table,2, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6','Taxonomy7','Taxonomy8','Taxonomy9'), sep = ';')

# find if there are non-euks here
i <- unique(otu.table$Taxonomy1)
# there is only "No blast hit" Remove these
j <- which(grepl(i[2],otu.table$Taxonomy1)) # there are 1683 rows with "No Blast Hit"
otu.table <- otu.table[-c(j),] 
# check
j <- which(grepl(i[2],otu.table$Taxonomy1)) # j is empty so it worked

#=====================================#
# Create Relative Abundance Matrices  #
#=====================================#

# First calculate relative abundance
otu.table.RA <- otu.table[,-c(2:12)]
otu.table.RA <- otu.table.RA[,-c(49)]
samplesums <- colSums(otu.table.RA[,-c(1)],na.rm = TRUE)
rm(i)

# Fill in RAs using loop
for(i in 2:ncol(otu.table.RA)){
  otu.table.RA[,i] <- otu.table.RA[,i]/samplesums[i-1]
}

# Remove low abundance/ bad sequencing samples
rm(i,j,k)
j <- unique(meta$'Eliminate?')
k <- which(grepl(j[2],meta$`Eliminate?`))
i <- meta$MariaKey[k]
i <- i[2:7] # NOTE- AE3a198A (i[1]) was already removed by Maria for not working
# So the abundance table doesn't have it even though it is in the meta data table
# Don't have to remove that one but remove others with low sequencing effort

# Need to convert colnames to number indexes in order to delete them
w <- which(colnames(otu.table.RA) %in% i)
otu.table.RA <- otu.table.RA[,-c(w)] 

# Check that all column sums are 1.0. For example:
colSums(otu.table.RA[,2:42], na.rm = TRUE)
# OK

# Calculate mean rel abun among sample duplicates
# Organization of samples is as follows:
# CAR212_PA = [AE3a103A AE3b103A AE3a198A AE3b198A AE3a234A AE3b234A AE3a295A AE3b295A AE3a314A AE3b314A AE3a900AM AE1b900AM];
# CAR212_FL = [AE3a103B AE3b103B AE3a198B AE3b198B AE3a234B AE3b234B AE3a295B AE3b295B AE3a314B AE3b314B AE3a900BM AE1b900BM];
# CAR216_PA = [AE2a143A AE2b143A AE2a200A AE2b200A AE2a237A AE2b237A AE2a247A AE2b247A AE2a267A AE2b267A AE2a900AN AE2b900AN];
# CAR216_FL = [AE2a143B AE2b143B AE2a200B AE2b200B AE2a237B AE2b237B AE2a247B AE2b247B AE2a267B AE2b267B AE2a900BN AE2b900BN];

# CAR212 is May. CAR216 is Nov. Capital A is PA size fraction. Capital B is FL.
# Depth is indicated by 3-digit number. Organized in order from shallowest 
# to deepest. 3 shallowest depths represent OXYCLINE  (O1, O2, O3). Next 2 depths are ANOXIC
# (A1, A2). Deepest depth is EUXINIC (E)

# Create matrix with means in columns in correct sample order 
C212A <- tibble("#OTU ID"=otu.table.RA$"#OTU ID", 
                    "O1"=rowMeans(otu.table.RA[ ,c("AE3a103A","AE3b103A")]),
                    "O2"=otu.table.RA$AE3b198A, # Sample AE3a198A was removed
                    "O3"=rowMeans(otu.table.RA[ ,c("AE3a234A","AE3b234A")]),
                    "A1"=rowMeans(otu.table.RA[ ,c("AE3a295A","AE3b295A")]),
                    "A2"=otu.table.RA$AE3a314A, # Sample AE3b314A was removed
                    "E"=rowMeans(otu.table.RA[ ,c("AE3a900AM","AE1b900AM")]))
C212B <- tibble("#OTU ID"=otu.table.RA$"#OTU ID", 
                    "O1"=rowMeans(otu.table.RA[ ,c("AE3a103B","AE3b103B")]),
                    "O2"=rowMeans(otu.table.RA[ ,c("AE3a198B","AE3b198B")]),
                    "O3"=rowMeans(otu.table.RA[ ,c("AE3a234B","AE3b234B")]),
                    "A1"=rowMeans(otu.table.RA[ ,c("AE3a295B","AE3b295B")]),
                    "A2"=rowMeans(otu.table.RA[ ,c("AE3a314B","AE3b314B")]),
                    "E"=rowMeans(otu.table.RA[ ,c("AE3a900BM","AE1b900BM")]))
C216A <- tibble("#OTU ID"=otu.table.RA$"#OTU ID", 
                    "O1"=rowMeans(otu.table.RA[ ,c("AE2a143A","AE2b143A")]),
                    "O2"=otu.table.RA$AE2b200A, # AE2a200A was removed
                    "O3"=rowMeans(otu.table.RA[ ,c("AE2a237A","AE2b237A")]),
                    "A1"=rowMeans(otu.table.RA[ ,c("AE2a247A","AE2b247A")]),
                    "A2"=rowMeans(otu.table.RA[ ,c("AE2a267A","AE2b267A")]),
                    "E"=otu.table.RA$AE2a900AN) # AE2b900AN was removed
C216B <- tibble("#OTU ID"=otu.table.RA$"#OTU ID", 
                    "O1"=rowMeans(otu.table.RA[ ,c("AE2a143B","AE2b143B")]),
                    "O2"=otu.table.RA$AE2b200B, # AE2a200B was removed
                    "O3"=rowMeans(otu.table.RA[ ,c("AE2a237B","AE2b237B")]),
                    "A1"=rowMeans(otu.table.RA[ ,c("AE2a247B","AE2b247B")]),
                    "A2"=otu.table.RA$AE2b267B, # AE2a267B was removed
                    "E"=otu.table.RA$AE2b900BN) # AE2a900BN was removed

# Check col sums. Should all be 1
colSums(C212A[, 2:7], na.rm = TRUE) # good
colSums(C212B[, 2:7], na.rm = TRUE) # good
colSums(C216A[, 2:7], na.rm = TRUE) # good
colSums(C216B[, 2:7], na.rm = TRUE) # good


# Re-link these relative abundace matrices with the taxanomic IDs in the same tibble
C212A <- left_join(C212A, otu.table[ , 1:10], by = "#OTU ID") 

C212B <- left_join(C212B, otu.table[ , 1:10], by = "#OTU ID") 

C216A <- left_join(C216A, otu.table[ , 1:10], by = "#OTU ID") 

C216B <- left_join(C216B, otu.table[ , 1:10], by = "#OTU ID") 


#====================================#
#====================================#
# Start with abundance bubble plots  #
#====================================#
#====================================#

######  Plot by Phylum ########
# High taxonomic level (~I am generally calling this phyla because it's one level below domain but it's really a mixture of
# phyla, superphyla, subgroups, because the taxonomy is so messed up even at this level)
tax2 <- unique(otu.table$Taxonomy2)
tax2 <- tax2[-c(10)] # removed level of uncategorized stuff at this level
tax2 <- sort(tax2)


# Make matrix for extracting each unique entry in tax2
O1_sums <- matrix(, nrow = length(tax2), ncol = 4)
O2_sums <- matrix(, nrow = length(tax2), ncol = 4)
O3_sums <- matrix(, nrow = length(tax2), ncol = 4)
A1_sums <- matrix(, nrow = length(tax2), ncol = 4)
A2_sums <- matrix(, nrow = length(tax2), ncol = 4)
E_sums <- matrix(, nrow = length(tax2), ncol = 4)


# Fill empty arrays with sums of that respective taxa and sample (O1, O2, etc)
for(i in 1:length(tax2)){
    x <- which(otu.table$Taxonomy2 == tax2[i]) # x is index (rows) of all instances of that taxon
      O1_sums[i,1] <- sum(C212A$O1[x], na.rm = TRUE)
      O2_sums[i,1] <- sum(C212A$O2[x], na.rm = TRUE)
      O3_sums[i,1] <- sum(C212A$O3[x], na.rm = TRUE)
      A1_sums[i,1] <- sum(C212A$A1[x], na.rm = TRUE)
      A2_sums[i,1] <- sum(C212A$A2[x], na.rm = TRUE)
      E_sums[i,1] <- sum(C212A$E[x], na.rm = TRUE)
      
      O1_sums[i,2] <- sum(C212B$O1[x], na.rm = TRUE)
      O2_sums[i,2] <- sum(C212B$O2[x], na.rm = TRUE)
      O3_sums[i,2] <- sum(C212B$O3[x], na.rm = TRUE)
      A1_sums[i,2] <- sum(C212B$A1[x], na.rm = TRUE)
      A2_sums[i,2] <- sum(C212B$A2[x], na.rm = TRUE)
      E_sums[i,2] <- sum(C212B$E[x], na.rm = TRUE)
      
      O1_sums[i,3] <- sum(C216A$O1[x], na.rm = TRUE)
      O2_sums[i,3] <- sum(C216A$O2[x], na.rm = TRUE)
      O3_sums[i,3] <- sum(C216A$O3[x], na.rm = TRUE)
      A1_sums[i,3] <- sum(C216A$A1[x], na.rm = TRUE)
      A2_sums[i,3] <- sum(C216A$A2[x], na.rm = TRUE)
      E_sums[i,3] <- sum(C216A$E[x], na.rm = TRUE)
      
      O1_sums[i,4] <- sum(C216B$O1[x], na.rm = TRUE)
      O2_sums[i,4] <- sum(C216B$O2[x], na.rm = TRUE)
      O3_sums[i,4] <- sum(C216B$O3[x], na.rm = TRUE)
      A1_sums[i,4] <- sum(C216B$A1[x], na.rm = TRUE)
      A2_sums[i,4] <- sum(C216B$A2[x], na.rm = TRUE)
      E_sums[i,4] <- sum(C216B$E[x], na.rm = TRUE)
      rm(x)
}

# Fill in tibbles
# Need a vector of OxyCond names that is 54 rows long, 1st 9 are 1st depth (O1), 2nd 10 are next depth (O2), etc
OxyCond <- paste(c("O1", "O1","O1","O1","O1","O1","O1","O1","O1","O2","O2","O2","O2","O2","O2","O2","O2","O2","O3","O3","O3","O3","O3","O3","O3","O3","O3", "A1","A1","A1","A1","A1","A1","A1","A1","A1","A2","A2","A2","A2","A2","A2","A2","A2","A2","E", "E", "E","E", "E","E","E","E","E"))
# Create tibble and fill with relative abundance data (only with C212A for now)
phylum.otu.table <- tibble("Taxonomy2" = paste(replicate(6, tax2), sep = " ", collapse = NULL), # concatenates the tax2 string vector 6 times
                     "OxyCond" = OxyCond,
                     "C212A" = c(O1_sums[1:9,1],O2_sums[1:9,1],O3_sums[1:9,1],A1_sums[1:9,1],A2_sums[1:9,1],E_sums[1:9,1]),
                     "C212B" = c(O1_sums[1:9,2],O2_sums[1:9,2],O3_sums[1:9,2],A1_sums[1:9,2],A2_sums[1:9,2],E_sums[1:9,2]),
                     "C216A" = c(O1_sums[1:9,3],O2_sums[1:9,3],O3_sums[1:9,3],A1_sums[1:9,3],A2_sums[1:9,3],E_sums[1:9,3]),
                     "C216B" = c(O1_sums[1:9,4],O2_sums[1:9,4],O3_sums[1:9,4],A1_sums[1:9,4],A2_sums[1:9,4],E_sums[1:9,4])
                    )

# Plot with ggplot 
#C212A
p1 = ggplot(phylum.otu.table,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy2))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(legend.position="none")

#C212B
p2 = ggplot(phylum.otu.table,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy2))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216A
p3 = ggplot(phylum.otu.table,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy2))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216B
p4 = ggplot(phylum.otu.table,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy2))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(axis.text.y=element_blank())

# Put all 4 panels on one plot
phylumfigure <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(2.83/10, 1.92/10, 1.92/10, 3.33/10)) # Played around with the ratios until panels looked even
phylumfigure <- annotate_figure(phylumfigure,
                bottom = text_grob("Oxygen Condition", hjust = 0.5, x = .5, size = 12))

# ggsave(filename = "Figures/Euk_phyla_bubbleplot.jpg", plot = phylumfigure, units = c("in"), width = 12, height = 3, dpi = 300, )
# ggsave(filename = "Figures/Euk_phyla_bubbleplot.eps", plot = phylumfigure, units = c("in"), width = 12, height = 3, dpi = 300, )







######  Plot Alveolata ########
# Plot all alveolates at high taxonomic resolution

# First need to pull out all alveolata from original otu table
Opisthokonta_otu_table <- otu.table %>%
  filter(Taxonomy2 == "Alveolata")

# For now try plotting down to Tax6, which gives group and clade of Syndiniales
# First make an additional column in otu table to retain full taxonomy, only up to Tax6, for plotting later:
Opisthokonta_otu_table[,ncol(Opisthokonta_otu_table)+1] <- paste(Opisthokonta_otu_table$Taxonomy3, ";", Opisthokonta_otu_table$Taxonomy4, ";", Opisthokonta_otu_table$Taxonomy5, ";", Opisthokonta_otu_table$Taxonomy6)
colnames(Opisthokonta_otu_table)[length(Opisthokonta_otu_table)] <- "Full_Taxonomy"

# Unique taxa at Tax level:
tax6 <- unique(Opisthokonta_otu_table$Full_Taxonomy)
tax6 <- sort(tax6)
tax6 # This gives 150 unique taxa


# Make matrix for extracting each unique entry in tax6
O1_sums <- matrix(, nrow = length(tax6), ncol = 4)
O2_sums <- matrix(, nrow = length(tax6), ncol = 4)
O3_sums <- matrix(, nrow = length(tax6), ncol = 4)
A1_sums <- matrix(, nrow = length(tax6), ncol = 4)
A2_sums <- matrix(, nrow = length(tax6), ncol = 4)
E_sums <- matrix(, nrow = length(tax6), ncol = 4)


# Fill empty arrays with sums of that respective taxa and sample (O1, O2, etc)
for(i in 1:length(tax6)){
  x <- which(Opisthokonta_otu_table$Full_Taxonomy == tax6[i]) # x is index (rows) of all instances of that taxon
  O1_sums[i,1] <- sum(C212A$O1[x], na.rm = TRUE)
  O2_sums[i,1] <- sum(C212A$O2[x], na.rm = TRUE)
  O3_sums[i,1] <- sum(C212A$O3[x], na.rm = TRUE)
  A1_sums[i,1] <- sum(C212A$A1[x], na.rm = TRUE)
  A2_sums[i,1] <- sum(C212A$A2[x], na.rm = TRUE)
  E_sums[i,1] <- sum(C212A$E[x], na.rm = TRUE)
  
  O1_sums[i,2] <- sum(C212B$O1[x], na.rm = TRUE)
  O2_sums[i,2] <- sum(C212B$O2[x], na.rm = TRUE)
  O3_sums[i,2] <- sum(C212B$O3[x], na.rm = TRUE)
  A1_sums[i,2] <- sum(C212B$A1[x], na.rm = TRUE)
  A2_sums[i,2] <- sum(C212B$A2[x], na.rm = TRUE)
  E_sums[i,2] <- sum(C212B$E[x], na.rm = TRUE)
  
  O1_sums[i,3] <- sum(C216A$O1[x], na.rm = TRUE)
  O2_sums[i,3] <- sum(C216A$O2[x], na.rm = TRUE)
  O3_sums[i,3] <- sum(C216A$O3[x], na.rm = TRUE)
  A1_sums[i,3] <- sum(C216A$A1[x], na.rm = TRUE)
  A2_sums[i,3] <- sum(C216A$A2[x], na.rm = TRUE)
  E_sums[i,3] <- sum(C216A$E[x], na.rm = TRUE)
  
  O1_sums[i,4] <- sum(C216B$O1[x], na.rm = TRUE)
  O2_sums[i,4] <- sum(C216B$O2[x], na.rm = TRUE)
  O3_sums[i,4] <- sum(C216B$O3[x], na.rm = TRUE)
  A1_sums[i,4] <- sum(C216B$A1[x], na.rm = TRUE)
  A2_sums[i,4] <- sum(C216B$A2[x], na.rm = TRUE)
  E_sums[i,4] <- sum(C216B$E[x], na.rm = TRUE)
  rm(x)
}


# Fill in tibbles
# Need a vector of OxyCond names that is 900 rows long, 1st 150 are 1st depth (O1), 2nd 150 are next depth (O2), etc
OxyCond <- c(replicate(length(tax6), c("O1")), replicate(length(tax6), c("O2")), replicate(length(tax6), c("O3")), 
              replicate(length(tax6), c("A1")),replicate(length(tax6), c("A2")), replicate(length(tax6), c("E")))
# Create tibble and fill with relative abundance data 
alveolata.otu.table <- tibble("Taxonomy" = paste(replicate(6, tax6), sep = " ", collapse = NULL), # concatenates the tax6 string vector 6 times
                           "OxyCond" = OxyCond,
                           "C212A" = c(O1_sums[1:length(tax6),1],O2_sums[1:length(tax6),1],O3_sums[1:length(tax6),1],A1_sums[1:length(tax6),1],A2_sums[1:length(tax6),1],E_sums[1:length(tax6),1]),
                           "C212B" = c(O1_sums[1:length(tax6),2],O2_sums[1:length(tax6),2],O3_sums[1:length(tax6),2],A1_sums[1:length(tax6),2],A2_sums[1:length(tax6),2],E_sums[1:length(tax6),2]),
                           "C216A" = c(O1_sums[1:length(tax6),3],O2_sums[1:length(tax6),3],O3_sums[1:length(tax6),3],A1_sums[1:length(tax6),3],A2_sums[1:length(tax6),3],E_sums[1:length(tax6),3]),
                           "C216B" = c(O1_sums[1:length(tax6),4],O2_sums[1:length(tax6),4],O3_sums[1:length(tax6),4],A1_sums[1:length(tax6),4],A2_sums[1:length(tax6),4],E_sums[1:length(tax6),4])
)

# Filter out those with abundance <0.1% (fractional abundance < 0.001)
# Write loop
# If, for a given tax6 (1-150), and given OxyCond (O1, O2, O3, A1, A2, E), there is a value in columns C212A, 
# C212B, C216A, or C216B for which the value is >0.001, retain all instances of that tax6

alveolata.otu.table.001 <- alveolata.otu.table[0,]
i <- 1
for(i in 1:length(tax6)){
  temp <- alveolata.otu.table %>%
    filter(Taxonomy==tax6[i])
  if(any(temp[,3:6] >= 0.001)) alveolata.otu.table.001 <- full_join(alveolata.otu.table.001, temp)
}






# Plot with ggplot 
#C212A
p1 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(legend.position="none")

#C212B
p2 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216A
p3 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216B
p4 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))+
  theme(axis.text.y=element_blank())

# Put all 4 panels on one plot
alveolatafigure_tax6 <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(51.5/100, 12.5/100, 12.5/100, 23.5/100)) # played with relative widths again because taxa labels are so long

alveolatafigure_tax6 <- annotate_figure(alveolatafigure_tax6,
                                bottom = text_grob("Oxygen Condition", hjust = 0.5, x = .5, size = 12))

# ggsave(filename = "Figures/Alveolata_tax6_bubbleplot.jpg", plot = alveolatafigure_tax6, units = c("in"), width = 16, height = 12, dpi = 300, )
# ggsave(filename = "Figures/Alveolata_tax6_bubbleplot.eps", plot = alveolatafigure_tax6, units = c("in"), width = 16, height = 12, dpi = 300, )
# This last figure is good becuase it gives a lot of detail on the individual clades of syndiniales but its
# prob more efficient for paper/ presetation to group at tax5, which gives just the groups (not clades) of syndiniales







# Plot Alveolates down to tax5
# First need to pull out all alveolata from original otu table
Alveolata_otu_table <- otu.table %>%
  filter(Taxonomy2 == "Alveolata")

# First make an additional column in otu table to retain full taxonomy, up to Tax5:
Alveolata_otu_table[,ncol(Alveolata_otu_table)+1] <- paste(Alveolata_otu_table$Taxonomy3, ";", Alveolata_otu_table$Taxonomy4, ";", Alveolata_otu_table$Taxonomy5)
colnames(Alveolata_otu_table)[length(Alveolata_otu_table)] <- "Full_Taxonomy"

# Unique taxa at Tax level:
tax5 <- unique(Alveolata_otu_table$Full_Taxonomy)
tax5 <- sort(tax5)
tax5 # This gives 54


# Make matrix for extracting each unique entry in tax5
O1_sums <- matrix(, nrow = length(tax5), ncol = 4)
O2_sums <- matrix(, nrow = length(tax5), ncol = 4)
O3_sums <- matrix(, nrow = length(tax5), ncol = 4)
A1_sums <- matrix(, nrow = length(tax5), ncol = 4)
A2_sums <- matrix(, nrow = length(tax5), ncol = 4)
E_sums <- matrix(, nrow = length(tax5), ncol = 4)


# Fill empty arrays with sums of that respective taxa and sample (O1, O2, etc)
for(i in 1:length(tax5)){
  x <- which(Alveolata_otu_table$Full_Taxonomy == tax5[i]) # x is index (rows) of all instances of that taxon
  O1_sums[i,1] <- sum(C212A$O1[x], na.rm = TRUE)
  O2_sums[i,1] <- sum(C212A$O2[x], na.rm = TRUE)
  O3_sums[i,1] <- sum(C212A$O3[x], na.rm = TRUE)
  A1_sums[i,1] <- sum(C212A$A1[x], na.rm = TRUE)
  A2_sums[i,1] <- sum(C212A$A2[x], na.rm = TRUE)
  E_sums[i,1] <- sum(C212A$E[x], na.rm = TRUE)
  
  O1_sums[i,2] <- sum(C212B$O1[x], na.rm = TRUE)
  O2_sums[i,2] <- sum(C212B$O2[x], na.rm = TRUE)
  O3_sums[i,2] <- sum(C212B$O3[x], na.rm = TRUE)
  A1_sums[i,2] <- sum(C212B$A1[x], na.rm = TRUE)
  A2_sums[i,2] <- sum(C212B$A2[x], na.rm = TRUE)
  E_sums[i,2] <- sum(C212B$E[x], na.rm = TRUE)
  
  O1_sums[i,3] <- sum(C216A$O1[x], na.rm = TRUE)
  O2_sums[i,3] <- sum(C216A$O2[x], na.rm = TRUE)
  O3_sums[i,3] <- sum(C216A$O3[x], na.rm = TRUE)
  A1_sums[i,3] <- sum(C216A$A1[x], na.rm = TRUE)
  A2_sums[i,3] <- sum(C216A$A2[x], na.rm = TRUE)
  E_sums[i,3] <- sum(C216A$E[x], na.rm = TRUE)
  
  O1_sums[i,4] <- sum(C216B$O1[x], na.rm = TRUE)
  O2_sums[i,4] <- sum(C216B$O2[x], na.rm = TRUE)
  O3_sums[i,4] <- sum(C216B$O3[x], na.rm = TRUE)
  A1_sums[i,4] <- sum(C216B$A1[x], na.rm = TRUE)
  A2_sums[i,4] <- sum(C216B$A2[x], na.rm = TRUE)
  E_sums[i,4] <- sum(C216B$E[x], na.rm = TRUE)
  rm(x)
}


# Fill in tibbles
# Need a vector of OxyCond names that is 900 rows long, 1st 150 are 1st depth (O1), 2nd 150 are next depth (O2), etc
OxyCond <- c(replicate(length(tax5), c("O1")), replicate(length(tax5), c("O2")), replicate(length(tax5), c("O3")), 
             replicate(length(tax5), c("A1")),replicate(length(tax5), c("A2")), replicate(length(tax5), c("E")))
# Create tibble and fill with relative abundance data 
alveolata.otu.table <- tibble("Taxonomy" = paste(replicate(6, tax5), sep = " ", collapse = NULL), # concatenates the tax5 string vector 6 times
                              "OxyCond" = OxyCond,
                              "C212A" = c(O1_sums[1:length(tax5),1],O2_sums[1:length(tax5),1],O3_sums[1:length(tax5),1],A1_sums[1:length(tax5),1],A2_sums[1:length(tax5),1],E_sums[1:length(tax5),1]),
                              "C212B" = c(O1_sums[1:length(tax5),2],O2_sums[1:length(tax5),2],O3_sums[1:length(tax5),2],A1_sums[1:length(tax5),2],A2_sums[1:length(tax5),2],E_sums[1:length(tax5),2]),
                              "C216A" = c(O1_sums[1:length(tax5),3],O2_sums[1:length(tax5),3],O3_sums[1:length(tax5),3],A1_sums[1:length(tax5),3],A2_sums[1:length(tax5),3],E_sums[1:length(tax5),3]),
                              "C216B" = c(O1_sums[1:length(tax5),4],O2_sums[1:length(tax5),4],O3_sums[1:length(tax5),4],A1_sums[1:length(tax5),4],A2_sums[1:length(tax5),4],E_sums[1:length(tax5),4])
)

# Filter out those with abundance <0.1% (fractional abundance < 0.001)
alveolata.otu.table.001 <- alveolata.otu.table[0,]
for(i in 1:length(tax5)){
  temp <- alveolata.otu.table %>%
    filter(Taxonomy==tax5[i])
  if(any(temp[,3:6] >= 0.001)) alveolata.otu.table.001 <- full_join(alveolata.otu.table.001, temp)
}






# Plot with ggplot 
#C212A
p1 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")

#C212B
p2 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216A
p3 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216B
p4 = ggplot(alveolata.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.25,.5,.75,1), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size=10))+
  theme(axis.text.y=element_blank())

# Put all 4 panels on one plot
alveolatafigure_tax5 <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(37/100, 16.5/100, 16.5/100, 30/100)) # played with relative widths again because taxa labels are so long

alveolatafigure_tax5 <- annotate_figure(alveolatafigure_tax5,
                                        bottom = text_grob("Oxygen Condition", hjust = 0.5, x = .5, size = 10))

# ggsave(filename = "Figures/Alveolata_tax5_bubbleplot.jpg", plot = alveolatafigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )
# ggsave(filename = "Figures/Alveolata_tax5_bubbleplot.eps", plot = alveolatafigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )







######  Plot Opisthokonta ########
# Plot all Opisthokonta at high taxonomic resolution
# Plot Alveolates down to tax5
# First need to pull out all alveolata from original otu table
Opisthokonta_otu_table <- otu.table %>%
  filter(Taxonomy2 == "Opisthokonta")

# First make an additional column in otu table to retain full taxonomy, up to Tax5:
Opisthokonta_otu_table[,ncol(Opisthokonta_otu_table)+1] <- paste(Opisthokonta_otu_table$Taxonomy3, ";", Opisthokonta_otu_table$Taxonomy4, ";", Opisthokonta_otu_table$Taxonomy5)
colnames(Opisthokonta_otu_table)[length(Opisthokonta_otu_table)] <- "Full_Taxonomy"

# Unique taxa at Tax level:
tax5 <- unique(Opisthokonta_otu_table$Full_Taxonomy)
tax5 <- sort(tax5)
tax5 # This gives 55


# Make matrix for extracting each unique entry in tax5
O1_sums <- matrix(, nrow = length(tax5), ncol = 4)
O2_sums <- matrix(, nrow = length(tax5), ncol = 4)
O3_sums <- matrix(, nrow = length(tax5), ncol = 4)
A1_sums <- matrix(, nrow = length(tax5), ncol = 4)
A2_sums <- matrix(, nrow = length(tax5), ncol = 4)
E_sums <- matrix(, nrow = length(tax5), ncol = 4)


# Fill empty arrays with sums of that respective taxa and sample (O1, O2, etc)
for(i in 1:length(tax5)){
  x <- which(Opisthokonta_otu_table$Full_Taxonomy == tax5[i]) # x is index (rows) of all instances of that taxon
  O1_sums[i,1] <- sum(C212A$O1[x], na.rm = TRUE)
  O2_sums[i,1] <- sum(C212A$O2[x], na.rm = TRUE)
  O3_sums[i,1] <- sum(C212A$O3[x], na.rm = TRUE)
  A1_sums[i,1] <- sum(C212A$A1[x], na.rm = TRUE)
  A2_sums[i,1] <- sum(C212A$A2[x], na.rm = TRUE)
  E_sums[i,1] <- sum(C212A$E[x], na.rm = TRUE)
  
  O1_sums[i,2] <- sum(C212B$O1[x], na.rm = TRUE)
  O2_sums[i,2] <- sum(C212B$O2[x], na.rm = TRUE)
  O3_sums[i,2] <- sum(C212B$O3[x], na.rm = TRUE)
  A1_sums[i,2] <- sum(C212B$A1[x], na.rm = TRUE)
  A2_sums[i,2] <- sum(C212B$A2[x], na.rm = TRUE)
  E_sums[i,2] <- sum(C212B$E[x], na.rm = TRUE)
  
  O1_sums[i,3] <- sum(C216A$O1[x], na.rm = TRUE)
  O2_sums[i,3] <- sum(C216A$O2[x], na.rm = TRUE)
  O3_sums[i,3] <- sum(C216A$O3[x], na.rm = TRUE)
  A1_sums[i,3] <- sum(C216A$A1[x], na.rm = TRUE)
  A2_sums[i,3] <- sum(C216A$A2[x], na.rm = TRUE)
  E_sums[i,3] <- sum(C216A$E[x], na.rm = TRUE)
  
  O1_sums[i,4] <- sum(C216B$O1[x], na.rm = TRUE)
  O2_sums[i,4] <- sum(C216B$O2[x], na.rm = TRUE)
  O3_sums[i,4] <- sum(C216B$O3[x], na.rm = TRUE)
  A1_sums[i,4] <- sum(C216B$A1[x], na.rm = TRUE)
  A2_sums[i,4] <- sum(C216B$A2[x], na.rm = TRUE)
  E_sums[i,4] <- sum(C216B$E[x], na.rm = TRUE)
  rm(x)
}


# Fill in tibbles
# Need a vector of OxyCond names that is 900 rows long, 1st 150 are 1st depth (O1), 2nd 150 are next depth (O2), etc
OxyCond <- c(replicate(length(tax5), c("O1")), replicate(length(tax5), c("O2")), replicate(length(tax5), c("O3")), 
             replicate(length(tax5), c("A1")),replicate(length(tax5), c("A2")), replicate(length(tax5), c("E")))
# Create tibble and fill with relative abundance data 
opisthokonta.otu.table <- tibble("Taxonomy" = paste(replicate(6, tax5), sep = " ", collapse = NULL), # concatenates the tax5 string vector 6 times
                              "OxyCond" = OxyCond,
                              "C212A" = c(O1_sums[1:length(tax5),1],O2_sums[1:length(tax5),1],O3_sums[1:length(tax5),1],A1_sums[1:length(tax5),1],A2_sums[1:length(tax5),1],E_sums[1:length(tax5),1]),
                              "C212B" = c(O1_sums[1:length(tax5),2],O2_sums[1:length(tax5),2],O3_sums[1:length(tax5),2],A1_sums[1:length(tax5),2],A2_sums[1:length(tax5),2],E_sums[1:length(tax5),2]),
                              "C216A" = c(O1_sums[1:length(tax5),3],O2_sums[1:length(tax5),3],O3_sums[1:length(tax5),3],A1_sums[1:length(tax5),3],A2_sums[1:length(tax5),3],E_sums[1:length(tax5),3]),
                              "C216B" = c(O1_sums[1:length(tax5),4],O2_sums[1:length(tax5),4],O3_sums[1:length(tax5),4],A1_sums[1:length(tax5),4],A2_sums[1:length(tax5),4],E_sums[1:length(tax5),4])
)

# Filter out those with abundance <0.1% (fractional abundance < 0.001)
opisthokonta.otu.table.001 <- opisthokonta.otu.table[0,]
for(i in 1:length(tax5)){
  temp <- opisthokonta.otu.table %>%
    filter(Taxonomy==tax5[i])
  if(any(temp[,3:6] >= 0.001)) opisthokonta.otu.table.001 <- full_join(opisthokonta.otu.table.001, temp)
}






# Plot with ggplot 
#C212A
p1 = ggplot(opisthokonta.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")

#C212B
p2 = ggplot(opisthokonta.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216A
p3 = ggplot(opisthokonta.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216B
p4 = ggplot(opisthokonta.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size=10))+
  theme(axis.text.y=element_blank())

# Put all 4 panels on one plot
opithokontafigure_tax5 <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(37/100, 16.5/100, 16.5/100, 30/100)) # played with relative widths again because taxa labels are so long

opithokontafigure_tax5 <- annotate_figure(opithokontafigure_tax5,
                                        bottom = text_grob("Oxygen Condition", hjust = 0.5, x = .5, size = 10))

# ggsave(filename = "Figures/Opisthokonta_tax5_bubbleplot.jpg", plot = opithokontafigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )
# ggsave(filename = "Figures/Opisthokonta_tax5_bubbleplot.eps", plot = opithokontafigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )








######  Plot Rhizaria ########
# Plot all Rhizaria at high taxonomic resolution
# Plot Alveolates down to tax5
# First need to pull out all alveolata from original otu table
Rhizaria_otu_table <- otu.table %>%
  filter(Taxonomy2 == "Rhizaria")

# First make an additional column in otu table to retain full taxonomy, up to Tax5:
Rhizaria_otu_table[,ncol(Rhizaria_otu_table)+1] <- paste(Rhizaria_otu_table$Taxonomy3, ";", Rhizaria_otu_table$Taxonomy4, ";", Rhizaria_otu_table$Taxonomy5)
colnames(Rhizaria_otu_table)[length(Rhizaria_otu_table)] <- "Full_Taxonomy"

# Unique taxa at Tax level:
tax5 <- unique(Rhizaria_otu_table$Full_Taxonomy)
tax5 <- sort(tax5)
tax5 # This gives 55


# Make matrix for extracting each unique entry in tax5
O1_sums <- matrix(, nrow = length(tax5), ncol = 4)
O2_sums <- matrix(, nrow = length(tax5), ncol = 4)
O3_sums <- matrix(, nrow = length(tax5), ncol = 4)
A1_sums <- matrix(, nrow = length(tax5), ncol = 4)
A2_sums <- matrix(, nrow = length(tax5), ncol = 4)
E_sums <- matrix(, nrow = length(tax5), ncol = 4)


# Fill empty arrays with sums of that respective taxa and sample (O1, O2, etc)
for(i in 1:length(tax5)){
  x <- which(Rhizaria_otu_table$Full_Taxonomy == tax5[i]) # x is index (rows) of all instances of that taxon
  O1_sums[i,1] <- sum(C212A$O1[x], na.rm = TRUE)
  O2_sums[i,1] <- sum(C212A$O2[x], na.rm = TRUE)
  O3_sums[i,1] <- sum(C212A$O3[x], na.rm = TRUE)
  A1_sums[i,1] <- sum(C212A$A1[x], na.rm = TRUE)
  A2_sums[i,1] <- sum(C212A$A2[x], na.rm = TRUE)
  E_sums[i,1] <- sum(C212A$E[x], na.rm = TRUE)
  
  O1_sums[i,2] <- sum(C212B$O1[x], na.rm = TRUE)
  O2_sums[i,2] <- sum(C212B$O2[x], na.rm = TRUE)
  O3_sums[i,2] <- sum(C212B$O3[x], na.rm = TRUE)
  A1_sums[i,2] <- sum(C212B$A1[x], na.rm = TRUE)
  A2_sums[i,2] <- sum(C212B$A2[x], na.rm = TRUE)
  E_sums[i,2] <- sum(C212B$E[x], na.rm = TRUE)
  
  O1_sums[i,3] <- sum(C216A$O1[x], na.rm = TRUE)
  O2_sums[i,3] <- sum(C216A$O2[x], na.rm = TRUE)
  O3_sums[i,3] <- sum(C216A$O3[x], na.rm = TRUE)
  A1_sums[i,3] <- sum(C216A$A1[x], na.rm = TRUE)
  A2_sums[i,3] <- sum(C216A$A2[x], na.rm = TRUE)
  E_sums[i,3] <- sum(C216A$E[x], na.rm = TRUE)
  
  O1_sums[i,4] <- sum(C216B$O1[x], na.rm = TRUE)
  O2_sums[i,4] <- sum(C216B$O2[x], na.rm = TRUE)
  O3_sums[i,4] <- sum(C216B$O3[x], na.rm = TRUE)
  A1_sums[i,4] <- sum(C216B$A1[x], na.rm = TRUE)
  A2_sums[i,4] <- sum(C216B$A2[x], na.rm = TRUE)
  E_sums[i,4] <- sum(C216B$E[x], na.rm = TRUE)
  rm(x)
}


# Fill in tibbles
OxyCond <- c(replicate(length(tax5), c("O1")), replicate(length(tax5), c("O2")), replicate(length(tax5), c("O3")), 
             replicate(length(tax5), c("A1")),replicate(length(tax5), c("A2")), replicate(length(tax5), c("E")))
# Create tibble and fill with relative abundance data 
rhizaria.otu.table <- tibble("Taxonomy" = paste(replicate(6, tax5), sep = " ", collapse = NULL), # concatenates the tax5 string vector 6 times
                                 "OxyCond" = OxyCond,
                                 "C212A" = c(O1_sums[1:length(tax5),1],O2_sums[1:length(tax5),1],O3_sums[1:length(tax5),1],A1_sums[1:length(tax5),1],A2_sums[1:length(tax5),1],E_sums[1:length(tax5),1]),
                                 "C212B" = c(O1_sums[1:length(tax5),2],O2_sums[1:length(tax5),2],O3_sums[1:length(tax5),2],A1_sums[1:length(tax5),2],A2_sums[1:length(tax5),2],E_sums[1:length(tax5),2]),
                                 "C216A" = c(O1_sums[1:length(tax5),3],O2_sums[1:length(tax5),3],O3_sums[1:length(tax5),3],A1_sums[1:length(tax5),3],A2_sums[1:length(tax5),3],E_sums[1:length(tax5),3]),
                                 "C216B" = c(O1_sums[1:length(tax5),4],O2_sums[1:length(tax5),4],O3_sums[1:length(tax5),4],A1_sums[1:length(tax5),4],A2_sums[1:length(tax5),4],E_sums[1:length(tax5),4])
)

# Filter out those with abundance <0.1% (fractional abundance < 0.001)
rhizaria.otu.table.001 <- rhizaria.otu.table[0,]
for(i in 1:length(tax5)){
  temp <- rhizaria.otu.table %>%
    filter(Taxonomy==tax5[i])
  if(any(temp[,3:6] >= 0.001)) rhizaria.otu.table.001 <- full_join(rhizaria.otu.table.001, temp)
}






# Plot with ggplot 
#C212A
p1 = ggplot(rhizaria.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")

#C212B
p2 = ggplot(rhizaria.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216A
p3 = ggplot(rhizaria.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216B
p4 = ggplot(rhizaria.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size=10))+
  theme(axis.text.y=element_blank())

# Put all 4 panels on one plot
rhizariafigure_tax5 <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(40/100, 15.5/100, 15.5/100, 29/100)) # played with relative widths again because taxa labels are so long

rhizariafigure_tax5 <- annotate_figure(rhizariafigure_tax5,
                                          bottom = text_grob("Oxygen Condition", hjust = 0.5, x = .5, size = 10))

# ggsave(filename = "Figures/Rhizaria_tax5_bubbleplot.jpg", plot = rhizariafigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )
# ggsave(filename = "Figures/Rhizaria_tax5_bubbleplot.eps", plot = rhizariafigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )






######  Plot Stramenopiles ########
# Plot all Stramenopiles at high taxonomic resolution
# Plot Stramenopiles down to tax5
# First need to pull out all  from original otu table
Stramenopiles_otu_table <- otu.table %>%
  filter(Taxonomy2 == "Stramenopiles")

# First make an additional column in otu table to retain full taxonomy, up to Tax5:
Stramenopiles_otu_table[,ncol(Stramenopiles_otu_table)+1] <- paste(Stramenopiles_otu_table$Taxonomy3, ";", Stramenopiles_otu_table$Taxonomy4, ";", Stramenopiles_otu_table$Taxonomy5)
colnames(Stramenopiles_otu_table)[length(Stramenopiles_otu_table)] <- "Full_Taxonomy"

# Unique taxa at Tax level:
tax5 <- unique(Stramenopiles_otu_table$Full_Taxonomy)
tax5 <- sort(tax5)
tax5 # This gives 58


# Make matrix for extracting each unique entry in tax5
O1_sums <- matrix(, nrow = length(tax5), ncol = 4)
O2_sums <- matrix(, nrow = length(tax5), ncol = 4)
O3_sums <- matrix(, nrow = length(tax5), ncol = 4)
A1_sums <- matrix(, nrow = length(tax5), ncol = 4)
A2_sums <- matrix(, nrow = length(tax5), ncol = 4)
E_sums <- matrix(, nrow = length(tax5), ncol = 4)


# Fill empty arrays with sums of that respective taxa and sample (O1, O2, etc)
for(i in 1:length(tax5)){
  x <- which(Stramenopiles_otu_table$Full_Taxonomy == tax5[i]) # x is index (rows) of all instances of that taxon
  O1_sums[i,1] <- sum(C212A$O1[x], na.rm = TRUE)
  O2_sums[i,1] <- sum(C212A$O2[x], na.rm = TRUE)
  O3_sums[i,1] <- sum(C212A$O3[x], na.rm = TRUE)
  A1_sums[i,1] <- sum(C212A$A1[x], na.rm = TRUE)
  A2_sums[i,1] <- sum(C212A$A2[x], na.rm = TRUE)
  E_sums[i,1] <- sum(C212A$E[x], na.rm = TRUE)
  
  O1_sums[i,2] <- sum(C212B$O1[x], na.rm = TRUE)
  O2_sums[i,2] <- sum(C212B$O2[x], na.rm = TRUE)
  O3_sums[i,2] <- sum(C212B$O3[x], na.rm = TRUE)
  A1_sums[i,2] <- sum(C212B$A1[x], na.rm = TRUE)
  A2_sums[i,2] <- sum(C212B$A2[x], na.rm = TRUE)
  E_sums[i,2] <- sum(C212B$E[x], na.rm = TRUE)
  
  O1_sums[i,3] <- sum(C216A$O1[x], na.rm = TRUE)
  O2_sums[i,3] <- sum(C216A$O2[x], na.rm = TRUE)
  O3_sums[i,3] <- sum(C216A$O3[x], na.rm = TRUE)
  A1_sums[i,3] <- sum(C216A$A1[x], na.rm = TRUE)
  A2_sums[i,3] <- sum(C216A$A2[x], na.rm = TRUE)
  E_sums[i,3] <- sum(C216A$E[x], na.rm = TRUE)
  
  O1_sums[i,4] <- sum(C216B$O1[x], na.rm = TRUE)
  O2_sums[i,4] <- sum(C216B$O2[x], na.rm = TRUE)
  O3_sums[i,4] <- sum(C216B$O3[x], na.rm = TRUE)
  A1_sums[i,4] <- sum(C216B$A1[x], na.rm = TRUE)
  A2_sums[i,4] <- sum(C216B$A2[x], na.rm = TRUE)
  E_sums[i,4] <- sum(C216B$E[x], na.rm = TRUE)
  rm(x)
}


# Fill in tibbles
OxyCond <- c(replicate(length(tax5), c("O1")), replicate(length(tax5), c("O2")), replicate(length(tax5), c("O3")), 
             replicate(length(tax5), c("A1")),replicate(length(tax5), c("A2")), replicate(length(tax5), c("E")))
# Create tibble and fill with relative abundance data 
stramenopiles.otu.table <- tibble("Taxonomy" = paste(replicate(6, tax5), sep = " ", collapse = NULL), # concatenates the tax5 string vector 6 times
                             "OxyCond" = OxyCond,
                             "C212A" = c(O1_sums[1:length(tax5),1],O2_sums[1:length(tax5),1],O3_sums[1:length(tax5),1],A1_sums[1:length(tax5),1],A2_sums[1:length(tax5),1],E_sums[1:length(tax5),1]),
                             "C212B" = c(O1_sums[1:length(tax5),2],O2_sums[1:length(tax5),2],O3_sums[1:length(tax5),2],A1_sums[1:length(tax5),2],A2_sums[1:length(tax5),2],E_sums[1:length(tax5),2]),
                             "C216A" = c(O1_sums[1:length(tax5),3],O2_sums[1:length(tax5),3],O3_sums[1:length(tax5),3],A1_sums[1:length(tax5),3],A2_sums[1:length(tax5),3],E_sums[1:length(tax5),3]),
                             "C216B" = c(O1_sums[1:length(tax5),4],O2_sums[1:length(tax5),4],O3_sums[1:length(tax5),4],A1_sums[1:length(tax5),4],A2_sums[1:length(tax5),4],E_sums[1:length(tax5),4])
)

# Filter out those with abundance <0.1% (fractional abundance < 0.001)
stramenopiles.otu.table.001 <- stramenopiles.otu.table[0,]
for(i in 1:length(tax5)){
  temp <- stramenopiles.otu.table %>%
    filter(Taxonomy==tax5[i])
  if(any(temp[,3:6] >= 0.001)) stramenopiles.otu.table.001 <- full_join(stramenopiles.otu.table.001, temp)
}






# Plot with ggplot 
#C212A
p1 = ggplot(stramenopiles.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")

#C212B
p2 = ggplot(stramenopiles.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C212B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C212B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216A
p3 = ggplot(stramenopiles.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216A))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216A")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme(axis.text.y=element_blank())

#C216B
p4 = ggplot(stramenopiles.otu.table.001,aes (x = OxyCond, y = fct_rev(as_factor(Taxonomy))))+ # the fancy stuff around y (Taxonomy2) just helps to present it in reverse order in the plot (from top to btm alphabetically)
  geom_point(aes(size = C216B))+
  scale_x_discrete(limits=c("O1","O2","O3","A1","A2","E"))+
  scale_size_area(breaks = c(0,.005,.01,.015), max_size = 6)+
  xlab("")+
  ylab("")+
  labs(size="Relative Abundance")+
  ggtitle("C216B")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        plot.title = element_text(size=10))+
  theme(axis.text.y=element_blank())

# Put all 4 panels on one plot
stramenopilesfigure_tax5 <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(40/100, 15.5/100, 15.5/100, 29/100)) # played with relative widths again because taxa labels are so long

stramenopilesfigure_tax5 <- annotate_figure(stramenopilesfigure_tax5,
                                       bottom = text_grob("Oxygen Condition", hjust = 0.5, x = .5, size = 10))

# ggsave(filename = "Figures/Stramenopiles_tax5_bubbleplot.jpg", plot = stramenopilesfigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )
# ggsave(filename = "Figures/Stramenopiles_tax5_bubbleplot.eps", plot = stramenopilesfigure_tax5, units = c("in"), width = 12, height = 4, dpi = 300, )





#====================================#
#====================================#
######  Network Analysis ########
#====================================#
#====================================#

# # Make new otu table with all euk OTUs that are >= 0.0001 (greater than 0.01%) in at least one sample
# 
# # First look at sample names in order to set up logical:
# samplenames <-  colnames(otu.table.RA)
# samplenames <- samplenames[2:length(samplenames)]
# 
# # Set up logical to retain only those rows where at least one samples shows a number greater than/ equal to 0.0001
# otu.table.RA.0001 = filter(otu.table.RA, AE3b103B >= 0.0001| AE3b103A >= 0.0001| AE3b234A >= 0.0001|  AE3a295A >= 0.0001| AE2b247A >= 0.0001|  AE2a247A >= 0.0001| AE3b234B >= 0.0001| AE2a143B >= 0.0001| AE3a198B >= 0.0001| AE3b198B >= 0.0001| AE3a234B >= 0.0001| AE1b900BM >= 0.0001| AE3a103A >= 0.0001| AE3a234A >= 0.0001| AE2b900BN >= 0.0001| AE2b143A >= 0.0001|  AE3a103B >= 0.0001| AE2a237B >= 0.0001| AE3b198A >= 0.0001| AE2b200B >= 0.0001| AE2b200A >= 0.0001| AE2b143B >= 0.0001| AE1b900AM >= 0.0001| AE3b295A >= 0.0001| AE2b237B >= 0.0001| AE3a314A >= 0.0001| AE2b247B >= 0.0001|  AE3a314B >= 0.0001| AE2a143A >= 0.0001| AE2a237A >= 0.0001| AE3a900BM >= 0.0001| AE2b237A >= 0.0001| AE2a267A >= 0.0001| AE3a900AM >= 0.0001| AE3b314B >= 0.0001| AE2a247B >= 0.0001| AE2b267B >= 0.0001| AE2b267A >= 0.0001| AE2a900AN >= 0.0001| AE3b295B >= 0.0001| AE3a295B >= 0.0001)
# # gives 5466 unique OTUs. This is enough to build a network and check it out.
# 
# # First reattach full OTU names from PR2 annotation
# otu.table.RA.0001 <- inner_join(otu.table.RA.0001, PR2_annotation, by = "#OTU ID") 
# 
# # Based on advice from Ashley Shade's lab (https://github.com/ShadeLab/DataTimeNotes/blob/master/20150401_DataTime_CompareNetworkAnalysis.md), 
# # SparCC does a decent job 
# # and since I have used it before I will stick with it
# 
# # subset tibble to make a matrix with numeric only
# euks.0001 <- as.data.frame(otu.table.RA.0001)
# euks.0001 <- euks.0001[,-c(43, 44, 45)]
# 
# # Make text file to run SparCC on this .txt files in python 
# colnames(euks.0001)[1] <- "OTU_ID"
# write.table(euks.0001, "Network_PR2_OTUs/euks.0001.txt", sep = "\t",row.names = FALSE, quote=FALSE) # STOPPED HERE FO RPYTHON  1/21
# 
# 
######  Run SparCC in R? ########
### Momentarily I thought I would do this in R rather than python. Ended up taking too long to figure out but here are my notes from that ###
#  Download the R SparCC package:
# https://rdrr.io/github/MPBA/r-sparcc/
#library(devtools)
## install_github("MPBA/r-sparcc")
#library(rsparcc)
## Run SparCC
# euks.0001_sparcc <- sparcc(euks.0001[,2:42], max.iter = 20, th = 0.1, exiter = 10)
# euks.0001_sparcc$CORR
## This is actually comparing sample to sample I think

##### Implement SparCC in Python 
## See lab notebook for exact coding

## CONTINUE HERE WHEN SPARCC IS DONE RUNNING ON EUKARYOTIC MATRIX

# # Bring in the correlation and p-value matrices
# euks.0001.corr <- as.matrix(read.table("Network_PR2_OTUs/euks.0001_basis_corr/cor_sparcc.out", header = T, row.names = 1, sep = "\t"))
# euks.0001.oneside <- as.matrix(read.table("Network_PR2_OTUs/euks.0001_pvals/pvals_one_sided.txt", header = T, row.names = 1, sep = "\t"))
# 
# # Remove everything above the diagnol in each matrix so results are not repeated later on when these are linearized
# euks.0001.corr[upper.tri(euks.0001.corr)] <- NA
# euks.0001.oneside[upper.tri(euks.0001.oneside)] <- NA
# 
# 
# 
######  Visualize Correlation Matrix ########
# NOTE- this takes forever to render on large matrices. OK to skip
# #install.packages('corrplot')
# library(corrplot)
# 
# # Combine the correlation and p-value matrices
# euks.0001.corr.pval <- list(euks.0001.corr, euks.0001.oneside)
# 
# # Plot the correlogram showing only significant correlations
# corrplot(euks.0001.corr.pval[[1]], type = "lower", order = "FPC", p.mat = euks.0001.corr.pval[[2]], sig.level = 0.05, insig = "blank", tl.pos = "n", diag = F, cl.pos = "r")
# # think this worked but plot is very large- taking a long time to render.




######  Plotting Networks ########
# Followed tutorial from http://www.kateto.net/wp-content/uploads/2015/06/Polnet%202015%20Network%20Viz%20Tutorial%20-%20Ognyanova.pdf
# (Downloaded pdf, put it in Mendeley)
# Reference Ognyanova 2015

####
## Side note-
## there is a different package here that may be easier to use
## from https://stackoverflow.com/questions/12373821/correlation-matrix-to-build-networks
##install.packages("qgraph")
#library("qgraph")
####



# install.packages("igraph") 
# install.packages("network") 
# install.packages("sna") 
# install.packages("ndtv")
# install.packages("RColorBrewer")
# install.packages("extrafont") 


# Note: I was having trouble with installations and dependencies. 
# From reading online it may be because I recently updated my Mac OS to Mojave
# Apparently when you do that, the update deletes command line tools
# https://community.rstudio.com/t/r-and-r-studio-on-mojave/15891
# Also I needed to update my R version from 3.3 to 3.5




# ####
# ## Try Tutorial Dataset from Polnet 2015 to see format:
# # Tutorial data (linkedin pdf) is in Polnet2015 folder:
# nodes <- read.csv("Polnet2015/Data/Dataset1-Media-Example-NODES.csv", header=T, as.is=T) 
# links <- read.csv("Polnet2015/Data/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
# # Examine tutorial data
# head(nodes) 
# head(links) 
# nrow(nodes); length(unique(nodes$id)) 
# nrow(links); nrow(unique(links[,c("from", "to")]))
# # Make network object
# net <- graph.data.frame(links, nodes, directed=T)
# net
# plot(net)
# ####


###### Plot Cariaco Euk Network #####
# 
# # Make a "links" matrix (the edges)
# # This contains the correlation weight of all pairs of nodes
# # So I have to simplify my correlation matrix into a list of pairs
# euk.links <- melt(euks.0001.corr)
# colnames(euk.links)[3] <- "weight"
# colnames(euk.links)
# # Put pvalues in this list for now
# euk.pvals <- melt(euks.0001.oneside)
# euk.links <- inner_join(euk.links, euk.pvals)
# colnames(euk.links)[1] <- "from"
# colnames(euk.links)[2] <- "to"
# colnames(euk.links)[4] <- "pval_1side"
# colnames(euk.links)
# 
# # Filter out non signficant relationships. then remove "pval" column
# euk.links.sig <- filter(euk.links, pval_1side <=0.01)
# euk.links.sig <- euk.links.sig[-c(4)]
# 
# # Add in "type" column to show positive or negative relationship. Then just retain the absolute value of the correlation
# euk.links.sig$type <- ifelse(euk.links.sig$weight > 0, 'positive',
#                                 ifelse(euk.links.sig$weight < 0, 'negative',"NA")
#                                 )
# euk.links.sig$weight <- abs(euk.links.sig$weight)
# 
# # Remove those with correl coef less than 0.15
# euk.links.sig <- filter(euk.links.sig, weight >=0.15)
# # 394 edges retained
# 
# # Redefine the vertices (nodes) by just what's left in the matrix after trimming out non sig relationships
# x <- melt(t(euk.links.sig[1:2]))
# x <- data.frame(unique(x$value))
# count(x) # there are 314 unique OTU nodes left
# colnames(x) <- colnames(PR2_annotation[1])
# 
# # Join the retained OTU IDs to taxonomy- this is the new "nodes" dataset
# euk.nodes <- inner_join(x, PR2_annotation)
# euk.nodes <- euk_nodes[,1:2]
# euk.nodes[,3:11] <- separate(euk.nodes,2, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6','Taxonomy7','Taxonomy8','Taxonomy9'), sep = ';')
# euk.nodes <- euk.nodes[,-c(3:4)] # remove redundant columns
# 
# # For plotting later, I also need to associate a number to each unique taxon within each Taxonomy colum (do Taxonomy2-4 for now)
# # First Taxonomy_2 (super phylum)
# x <- 1:as.numeric(count(data.frame(unique(euk.nodes$Taxonomy2)))[1]) # there are 7 categories
# x_labs <- as.array(unique(euk.nodes$Taxonomy2))
# # Make new column for number index ot Taxonomy2 called Taxonomy2_label
# euk.nodes[,dim(euk.nodes)[2]+1] <- NA
# colnames(euk.nodes)[dim(euk.nodes)[2]] <- "Taxonomy2_label"
# # Fill in using loop
# for(j in 1:length(x)) {
# for(i in 1:dim(euk.nodes)[1]) {
#     if(euk.nodes$Taxonomy2[i] == x_labs[j]){euk.nodes$Taxonomy2_label[i] <- x[j]}
# }
# }
# 
# # Next Taxonomy_3 (phylum)
# rm(x, i, j)
# x <- 1:as.numeric(count(data.frame(unique(euk.nodes$Taxonomy3)))[1]) # there are 18 categories
# x_labs <- as.array(unique(euk.nodes$Taxonomy3))
# # Make new column for number index 
# euk.nodes[,dim(euk.nodes)[2]+1] <- NA
# colnames(euk.nodes)[dim(euk.nodes)[2]] <- "Taxonomy3_label"
# # Fill in using loop
# for(j in 1:length(x)) {
#   for(i in 1:dim(euk.nodes)[1]) {
#     if(euk.nodes$Taxonomy3[i] == x_labs[j]){euk.nodes$Taxonomy3_label[i] <- x[j]}
#   }
# }
# 
# 
# # Next Taxonomy_4 (class)
# rm(x, i, j)
# x <- 1:as.numeric(count(data.frame(unique(euk.nodes$Taxonomy4)))[1]) # there are 45 categories
# x_labs <- as.array(unique(euk.nodes$Taxonomy4))
# # Make new column for number index 
# euk.nodes[,dim(euk.nodes)[2]+1] <- NA
# colnames(euk.nodes)[dim(euk.nodes)[2]] <- "Taxonomy4_label"
# # Fill in using loop
# for(j in 1:length(x)) {
#   for(i in 1:dim(euk.nodes)[1]) {
#     if(euk.nodes$Taxonomy4[i] == x_labs[j]){euk.nodes$Taxonomy4_label[i] <- x[j]}
#   }
# }
# 
# # Next Taxonomy_5 (order?)
# rm(x, i, j)
# x <- 1:as.numeric(count(data.frame(unique(euk.nodes$Taxonomy5)))[1]) # there are 63 categories
# x_labs <- as.array(unique(euk.nodes$Taxonomy5))
# # Make new column for number index 
# euk.nodes[,dim(euk.nodes)[2]+1] <- NA
# colnames(euk.nodes)[dim(euk.nodes)[2]] <- "Taxonomy5_label"
# # Fill in using loop
# for(j in 1:length(x)) {
#   for(i in 1:dim(euk.nodes)[1]) {
#     if(euk.nodes$Taxonomy5[i] == x_labs[j]){euk.nodes$Taxonomy5_label[i] <- x[j]}
#   }
# }
# 
# # Next Taxonomy_6 (family?- this might be a mixed bag of levels)
# rm(x, i, j)
# x <- 1:as.numeric(count(data.frame(unique(euk.nodes$Taxonomy6)))[1]) # there are 103 categories
# x_labs <- as.array(unique(euk.nodes$Taxonomy6))
# # Make new column for number index 
# euk.nodes[,dim(euk.nodes)[2]+1] <- NA
# colnames(euk.nodes)[dim(euk.nodes)[2]] <- "Taxonomy6_label"
# # Fill in using loop
# for(j in 1:length(x)) {
#   for(i in 1:dim(euk.nodes)[1]) {
#     if(euk.nodes$Taxonomy6[i] == x_labs[j]){euk.nodes$Taxonomy6_label[i] <- x[j]}
#   }
# }
# 
# 
# 
# # Remake network object
# euk.net <- graph.data.frame(euk.links.sig, euk.nodes, directed=F)
# # check attributes of network object
# euk.net
# E(euk.net) # shows 197 edges
# V(euk.net) # there are 314 vertices
# E(euk.net)$type
# V(euk.net)$Taxonomy
# V(euk.net)$Taxonomy2
# V(euk.net)$Taxonomy2_label
# 
# # Plot
# plot(euk.net) 
# # Finally! Renders slowly but much faster than before
# 
# # Personalize plot
# # Generate colors based Taxonomy2. Use the brewer palletes
# unique(V(euk.net)$Taxonomy2) # there are 7 unique units
# colrs <- c(brewer.pal(7, "Dark2"))
# V(euk.net)$color <- colrs[V(euk.net)$Taxonomy2_label]
# 
# # Compute node degrees (#links) and use that to set node size: 
# deg <- igraph::degree(euk.net, mode="all") # had to specifiy igraph::degree because there is another function degree
# V(euk.net)$size <- deg^2 # multiplication factor to play with
# 
# # Turn off labels for now
# V(euk.net)$label <- NA
# 
# # Set edge width based on weight: 
# E(euk.net)$width <- E(euk.net)$weight*10
# 
# #change edge color: 
# E(euk.net)$edge.color <- "black" 
# #E(euk.net)$width <- 1+E(euk.net)$weight/12 
# 
# # plot
# plot(euk.net) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # Looks good
# 
# # Play around with layouts
# # random spread of nodes
# plot(euk.net, layout=layout.random) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # nodes in a circle
# plot(euk.net, layout=layout.circle(euk.net)) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # nodes randomly spread within a sphere
# plot(euk.net, layout=layout.sphere(euk.net)) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # Fruchterman-Reingold (a commonly used force-directed layout)
# # Some info about this layout from the tutorial: 
# # Force-directed layouts try to get a nice-looking graph where edges are similar in length and cross each other as little 
# # as possible. They simulate the graph as a physical system. Nodes are electrically charged particles that repulse each 
# # other when they get too close. The edges act as springs that attract connected nodes closer together. As a result, 
# # nodes are evenly distributed through the chart area, and the layout is intuitive in that nodes which share more 
# # connections are closer to each other. The disadvantage of these algorithms is that they are rather slow and therefore 
# # less often used in graphs larger than ~1000 vertices.
# # You can adjust some of the parameters of these- but i keep getting warnings when I do this
# # Also I think this is already the defauly layout
# l <- layout.fruchterman.reingold(euk.net, repulserad=vcount(euk.net)^.5, area=vcount(euk.net)^2)
# plot(euk.net, layout=l) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # Another layout is Kamada Kawai
# plot(euk.net, layout=layout.kamada.kawai(euk.net)) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # And Spring- but this is not working for me
# plot(euk.net, layout=layout.spring(euk.net, mass=.5)) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # And LGL, which is for highly connected networks
# plot(euk.net, layout=layout.lgl) 
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# # Something called Reingolf Tilford
# plot(euk.net, layout.reingold.tilford(euk.net, circular=T))
# legend(x=-1.5, y=-1.1, c(unique(V(euk.net)$Taxonomy2)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# 
# 
# # I don't love any of these but the best is probably the default Fruchterman-Reingold or the circular
# 
# # Write a loop to generate all layouts and pick best ones
# layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE) 
# # Remove layouts that do not apply to our graph. 
# layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]
# par(mfrow=c(3,3))
# for (layout in layouts) { 
#   print(layout) 
#   l <- do.call(layout, list(euk.net)) 
#   plot(euk.net, edge.arrow.mode=0, layout=l, main=layout)
#   }
# 
# dev.off()
# 
# # I exported these layouts to pdfs to take a look at them all together. See ppt
# # Also notice that "auto" picks the best one for the data, which is fruchterman.reingold
# # And I got warnings for layout.fruchterman.reingold.grid, 
# 
# # Upon looking at them all together, I think star or circle is the best. It minimizes the couple relationships and emphasizes the nodes that have multiple connections.
# 
# # As of 1/23 waiting for Sparcc to re-run network plots on euk dataset. SParcc files are being generated (in Amazon AWS).
# # Redid analysis with any OTU > 0.01% (instead of 0.1%). This created a huge dataset so SparCC is taking a while
# 
# 
# 
# 
# 








######  Make Master OTU Table and Run SparCC (3 domains) #####
# Make a master OTU table which includes relative abundances of all 3 domains, organized by sample (column) and OTU ID (row)
# as well as environmental parameters (metadata)
# I already have the euk table:
euk.otu.table <- otu.table.RA
euk.otu.table <- left_join(euk.otu.table, PR2_annotation, by = "#OTU ID") 
euk.otu.table <- euk.otu.table[-c(50:51)]
euk.otu.table <- separate(euk.otu.table,length(euk.otu.table), c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6','Taxonomy7','Taxonomy8','Taxonomy9', 'Taxonomy10'), 
                            sep = ";") # separate taxonomy
euk.otu.table[,49:58] <- lapply(euk.otu.table[,49:58], gsub, pattern='_', replacement=' ') # replace _ with spaces


# Read in raw abundance table with SILVA annotations
bac_abun_table <- read_delim("iTags_Analyses_SILVA_old/Cariaco_AB_updated_long_sums.txt",delim = "\t", col_names = T)
arch_abun_table <- read_delim("iTags_Analyses_SILVA_old/Cariaco_AA_updated_long_sums.txt",delim = "\t", col_names = T)


# separate taxonomy into multiple columns
# Bacteria
bac_abun_table <- separate(bac_abun_table,51, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6','Taxonomy7','Taxonomy8','Taxonomy9', 'Taxonomy10'), 
                           sep = ";")
# For some reason I can't separate by "__" so just remove double underscores from the names using lapply and gsub
bac_abun_table[,51:60] <- lapply(bac_abun_table[,51:60], gsub, pattern='__', replacement='')
# And replace single underscores with spaces
bac_abun_table[,51:60] <- lapply(bac_abun_table[,51:60], gsub, pattern='_', replacement=' ')
# Remove "sum" column
bac_abun_table <- bac_abun_table[-c(50)]


# Archaea
arch_abun_table <- separate(arch_abun_table,49, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6','Taxonomy7','Taxonomy8','Taxonomy9', 'Taxonomy10'), 
                           sep = ";")
arch_abun_table[,49:58] <- lapply(arch_abun_table[,49:58], gsub, pattern='__', replacement='')
# And replace single underscores with spaces
arch_abun_table[,49:58] <- lapply(arch_abun_table[,49:58], gsub, pattern='_', replacement=' ')
# Remove "sum" column
arch_abun_table <- arch_abun_table[-c(48)]


# There are also samples that we previously identified as bad due to low sequencing effort
# See "SampleKey.xlsx"
# Remove those
colSums(bac_abun_table[c("AB3a900A", "AB2a200A", "AB2b267A")],na.rm = TRUE)
bac_abun_table <- subset(bac_abun_table, select = -c(AB3a900A, AB2a200A, AB2b267A))

colSums(arch_abun_table[c("AA2b900AN", "AA2a247B", "AA2a900BN", "AA2b900BN")],na.rm = TRUE)
arch_abun_table <- subset(arch_abun_table, select = -c(AA2b900AN, AA2a247B, AA2a900BN, AA2b900BN))


# Calculate relative abundances 
# Bacteria
samplesums <- colSums(bac_abun_table[,2:46],na.rm = TRUE)
samplesums
bac.otu.table <- bac_abun_table #Use this as the empty matrix to fill in
bac.otu.table[,2:46] <- NA
# Fill in Relative abundances using loop
rm(i)
for(i in 2:46){
  bac.otu.table[,i] <- bac_abun_table[,i]/samplesums[i-1]
}
# Check column sums to make sure everything = 1
colSums(bac.otu.table[,2:46],na.rm = TRUE) # worked

# Archaea
rm(i, samplesums)
samplesums <- colSums(arch_abun_table[,2:43],na.rm = TRUE)
samplesums
arch.otu.table <- arch_abun_table #Use this as the empty matrix to fill in
arch.otu.table[,2:43] <- NA
# Fill in Relative abundances using loop
for(i in 2:43){
  arch.otu.table[,i] <- arch_abun_table[,i]/samplesums[i-1]
}
# Check column sums to make sure everything = 1
colSums(arch.otu.table[,2:43],na.rm = TRUE) # worked


# Put all Euks, Bac, and Arch in same table and concatenate by sample names

# First remove the prefix from sample names
# Eg Sample from 314m, fraction A for bacteria is called AB3a314A, for Archaea is called AA3a314A, and for Euks is called AE3a314A
colnames(euk.otu.table) <- lapply(colnames(euk.otu.table), gsub, pattern='AE', replacement='A')
colnames(bac.otu.table) <- lapply(colnames(bac.otu.table), gsub, pattern='AB', replacement='A')
colnames(arch.otu.table) <- lapply(colnames(arch.otu.table), gsub, pattern='AA', replacement='A')

# Make one correction to colnames
# the Deep (900m) samples from May in the bacteria and archaea dataset are called 3a900A, 1b900A, 3a900B, 1b900B
# while in the eukaryote dataset they are called 3a900AM, 1b900AM, 3a900BM, 1b900BM
# remove Ms from the names in euk dataset
euk.otu.table <- rename(euk.otu.table, "A3a900A" = "A3a900AM")
euk.otu.table <- rename(euk.otu.table, "A1b900A" = "A1b900AM")
euk.otu.table <- rename(euk.otu.table, "A3a900B" = "A3a900BM")
euk.otu.table <- rename(euk.otu.table, "A1b900B" = "A1b900BM")

colnames(euk.otu.table)

# Concatenate
temp <- full_join(arch.otu.table, bac.otu.table)
otu.table <- full_join(temp,euk.otu.table)
colnames(otu.table)

# Need to remove some rows. 
# Change name of "#OTU ID" because the # is messing up a lot of my commands
colnames(otu.table)[1] <- "OTU_ID"
colnames(euk.otu.table)[1] <- "OTU_ID"
colnames(bac.otu.table)[1] <- "OTU_ID"
colnames(arch.otu.table)[1] <- "OTU_ID"

# There are "sums" rows left over from some of the datasets and rows with "NA" across all columns 
filter(otu.table, !grepl('denovo', OTU_ID))
# Filter these non-data rows out
otu.table <- filter(otu.table, grepl('denovo', OTU_ID))

# Filter by most abundant. Retain any OTU that is at least 0.5% (RA = 0.005) in any sample
otu.table.005 <- filter(otu.table, A3a314A >= 0.005 | A2a237B >= 0.005 | A2b237A >= 0.005 | A3a900B >= 0.005 | A1b900A >= 0.005 | A3a103B >= 0.005 | 
                          A2a247A >= 0.005 | A1b900B >= 0.005 | A3b198A >= 0.005 | A2b237B >= 0.005 | A2b143A >= 0.005 | A3b103A >= 0.005 | 
                          A2a143A >= 0.005 | A3a103A>= 0.005 | A3b198B >= 0.005 | A3b234B >= 0.005 | A2a143B >= 0.005 | A3a234A >= 0.005 | 
                          A3a295B >= 0.005 | A2b247B >= 0.005 | A3a234B >= 0.005 | A2b247A >= 0.005 | A3b234A >= 0.005 | A3a314B >= 0.005 | 
                          A3a198B >= 0.005 | A3b295A >= 0.005 | A3a295A >= 0.005 | A2a237A >= 0.005 | A2b200A >= 0.005 | A2b267B >= 0.005 | 
                          A3b314B >= 0.005 | A2b143B >= 0.005 | A2b200B >= 0.005 | A3b103B >= 0.005 | A3b295B >= 0.005 | A2a900AN>= 0.005)



# Add in environmental data as rows
# upload metadata
meta <- read_delim("Metadata.txt",delim = "\t", col_names = T)

# Need to extract numeric metadata only
meta2 <- meta[6:23,]
meta2[,2:49] <- lapply(meta2[,2:49], as.numeric)

# Rename "key" to "OTU_ID" even though it's not really OTU- this will help match the dataframes
colnames(meta2)[1] <- colnames(otu.table.005)[1]

# Match to master OTU table abundant (>0.5%) OTUs
otu.table.005 <- full_join(otu.table.005, meta2)


# There are some libraries which worked for bacteria, for example, but not archaea. For network analysis need to remove 
# those samples which many not have corresponding values in all 3 datasets.
clean <- which(colnames(otu.table.005) %in% c("A3a900A", "A2a200A", "A2b267A", "A2a267A", "A2b900AN", "A2a247B", 
                                          "A2a900BN", "A2b900BN", "A3a198A", "A3b314A", "A2a200B", "A2a267B"))
otu.table.005 <- otu.table.005[, -clean] 




# There are 377 OTUs remaining.
which(otu.table.005$Taxonomy1 == "Archaea")
# Rows 1-60 are Archaea
which(otu.table.005$Taxonomy1 == "Bacteria")
# Rows 61-137 are Bacteria
which(otu.table.005$Taxonomy1 == "Eukaryota")
# Rows 138-377 are Eukaryota
# Meaning rows 378-395 are the metadata



# Export to txt file (leaving out strings except denovo ID)
write.table(otu.table.005[,1:37], "Network_PR2_OTUs/otu.table.005.txt", sep = "\t",row.names = FALSE, quote=FALSE)

# Run correlation analyses on these in SparCC in Python
# See lab notebook for code





##### Plot 3-Domain Network ####
# Adapting from Network code above for eukaryotes only

# Import results from SparCC
otu.table.005.corr <- as.matrix(read.table("Network_PR2_OTUs/otu.table.005_basis_corr/cor_sparcc.out", header = T, row.names = 1, sep = "\t"))
otu.table.005.oneside <- as.matrix(read.table("Network_PR2_OTUs/otu.table.005_pvals/pvals_one_sided.txt", header = T, row.names = 1, sep = "\t"))

# Remove everything above the diaganol in each matrix so results are not repeated later on when these are linearized
otu.table.005.corr[upper.tri(otu.table.005.corr)] <- NA
otu.table.005.oneside[upper.tri(otu.table.005.oneside)] <- NA

# Associate the denovo IDs with some marker that indicates which of the 3 domains 
#(because there is overlap of denove IDs from each dataset)
# I indicated above which of the 377 OTUs are from which domain (see L1658)
# Rows 1-60 are Archaea
# Rows 61-137 are Bacteria
# Rows 138-377 are Eukaryota
# Rows 378-395 are the metadata
colnames(otu.table.005.corr)[1:60] <- paste(colnames(otu.table.005.corr)[1:60],"-A")
rownames(otu.table.005.corr)[1:60] <- paste(rownames(otu.table.005.corr)[1:60],"-A")
colnames(otu.table.005.corr)[61:137] <- paste(colnames(otu.table.005.corr)[61:137],"-B")
rownames(otu.table.005.corr)[61:137] <- paste(rownames(otu.table.005.corr)[61:137],"-B")
colnames(otu.table.005.corr)[138:377] <- paste(colnames(otu.table.005.corr)[138:377],"-E")
rownames(otu.table.005.corr)[138:377] <- paste(rownames(otu.table.005.corr)[138:377],"-E")

colnames(otu.table.005.oneside)[1:60] <- paste(colnames(otu.table.005.oneside)[1:60],"-A")
rownames(otu.table.005.oneside)[1:60] <- paste(rownames(otu.table.005.oneside)[1:60],"-A")
colnames(otu.table.005.oneside)[61:137] <- paste(colnames(otu.table.005.oneside)[61:137],"-B")
rownames(otu.table.005.oneside)[61:137] <- paste(rownames(otu.table.005.oneside)[61:137],"-B")
colnames(otu.table.005.oneside)[138:377] <- paste(colnames(otu.table.005.oneside)[138:377],"-E")
rownames(otu.table.005.oneside)[138:377] <- paste(rownames(otu.table.005.oneside)[138:377],"-E")


# Make a "links" matrix (the edges) with correlation weights of all pairs of nodes
# Simplify correlation matrix into a list of pairs
otu.table.005.links <- melt(otu.table.005.corr)
colnames(otu.table.005.links)[3] <- "weight"
colnames(otu.table.005.links)
# Put pvalues in this list for now
otu.table.005.pvals <- melt(otu.table.005.oneside)
otu.table.005.links <- inner_join(otu.table.005.links, otu.table.005.pvals)
colnames(otu.table.005.links)[1] <- "from"
colnames(otu.table.005.links)[2] <- "to"
colnames(otu.table.005.links)[4] <- "pval_1side"
colnames(otu.table.005.links)


# Filter out non signficant relationships. then remove "pval" column
otu.table.005.links.sig <- filter(otu.table.005.links, pval_1side <=0.05)
otu.table.005.links.sig <- otu.table.005.links.sig[-c(4)]
# Reduced from 156,025 links to 9259


# Add in "type" column to show positive or negative relationship. Then just retain the absolute value of the correlation
otu.table.005.links.sig$type <- ifelse(otu.table.005.links.sig$weight > 0, 'positive',
                             ifelse(otu.table.005.links.sig$weight < 0, 'negative',"NA")
)
otu.table.005.links.sig$weight <- abs(otu.table.005.links.sig$weight)

# Remove those with correl coef less than 0.125
otu.table.005.links.sig <- filter(otu.table.005.links.sig, weight >=0.125)
# 577 edges retained
# The last relationships in the list are correlations between the environmental variables. I don't care about these for the figure. Remove them.
otu.table.005.links.sig <- otu.table.005.links.sig[-c(546:577), ]


# Redefine the vertices (nodes) by just what's left in the matrix after trimming out non sig relationships
x <- melt(t(otu.table.005.links.sig[1:2]))
x <- data.frame(unique(x$value))
count(x) # there are 168 unique OTUs still remaining
colnames(x) <- "OTU_ID"

# Join the retained OTU IDs to taxonomy- this is the new "nodes" dataset

# I think the easiest way to do this is take the taxonomy/ OTU ID from arch_abun_table, add in "-A" 
# to the denovo IDs, then match by left_join
arch_taxonomy <- arch.otu.table[-c(2:43)]
arch_taxonomy$OTU_ID <- paste(arch_taxonomy$OTU_ID,"-A")
bac_taxonomy <- bac.otu.table[-c(2:46)]
bac_taxonomy$OTU_ID <- paste(bac_taxonomy$OTU_ID,"-B")
euk_taxonomy <- euk.otu.table[-c(2:42)]
euk_taxonomy$OTU_ID <- paste(euk_taxonomy$OTU_ID,"-E")

taxonomy <- full_join(arch_taxonomy, bac_taxonomy) #721776 rows, 11 cols
taxonomy <- full_join(taxonomy, euk_taxonomy) # 1167818, 11 cols

otu.table.005.nodes <- left_join(x, taxonomy)

# For plotting later, I also need to associate a number to each unique taxon within each Taxonomy colum 
# Taxonomy_1 (domain)
rm(x,i,j)
x_labs <- as.array(unique(otu.table.005.nodes$Taxonomy1[!is.na(otu.table.005.nodes$Taxonomy1)])) # removing NA
x <- 1:as.numeric(count(data.frame(x_labs))) # there are 30 categories
# Make new column for number index of Taxonomy1 called Taxonomy1_label
otu.table.005.nodes[,dim(otu.table.005.nodes)[2]+1] <- NA
colnames(otu.table.005.nodes)[dim(otu.table.005.nodes)[2]] <- "Taxonomy1_label"
# Fill in using loop
for(j in 1:length(x)) {
  for(i in 1:dim(otu.table.005.nodes)[1]) {
    if(is.na(otu.table.005.nodes$Taxonomy1)[i]) next
    if(otu.table.005.nodes$Taxonomy1[i] == x_labs[j]){otu.table.005.nodes$Taxonomy1_label[i] <- x[j]}
  }
}


# Taxonomy_2 (phylum/ super phylum) 
# NOTE here there are NAs in taxononmy name so need to adjust accordingly otherwise loops don't work
rm(x,i,j)
x_labs <- as.array(unique(otu.table.005.nodes$Taxonomy2[!is.na(otu.table.005.nodes$Taxonomy2)])) # removing NA
x <- 1:as.numeric(count(data.frame(x_labs))) # there are 30 categories
# Make new column for number index of Taxonomy2 called Taxonomy2_label
otu.table.005.nodes[,dim(otu.table.005.nodes)[2]+1] <- NA
colnames(otu.table.005.nodes)[dim(otu.table.005.nodes)[2]] <- "Taxonomy2_label"
# Fill in using loop- had to modify this a bit for cases where there is "NA" in Taxonomy2
for(j in 1:length(x)) {
  for(i in 1:dim(otu.table.005.nodes)[1]) {
    if(is.na(otu.table.005.nodes$Taxonomy2)[i]) next
    if(otu.table.005.nodes$Taxonomy2[i] == x_labs[j]){otu.table.005.nodes$Taxonomy2_label[i] <- x[j]}
  }
}


# Next Taxonomy_3
rm(x,i,j)
x_labs <- as.array(unique(otu.table.005.nodes$Taxonomy3[!is.na(otu.table.005.nodes$Taxonomy3)])) # 
x <- 1:as.numeric(count(data.frame(x_labs))) # there are 72 categories
# Make new column for number index of Taxonomy2 called Taxonomy3_label
otu.table.005.nodes[,dim(otu.table.005.nodes)[2]+1] <- NA
colnames(otu.table.005.nodes)[dim(otu.table.005.nodes)[2]] <- "Taxonomy3_label"
# Fill in using loop- 
for(j in 1:length(x)) {
  for(i in 1:dim(otu.table.005.nodes)[1]) {
    if(is.na(otu.table.005.nodes$Taxonomy3)[i]) next
    if(otu.table.005.nodes$Taxonomy3[i] == x_labs[j]){otu.table.005.nodes$Taxonomy3_label[i] <- x[j]}
  }
}


# Next Taxonomy_4
rm(x,i,j)
x_labs <- as.array(unique(otu.table.005.nodes$Taxonomy4[!is.na(otu.table.005.nodes$Taxonomy4)])) # 
x <- 1:as.numeric(count(data.frame(x_labs))) # there are 138 categories
# Make new column for number index of Taxonomy2 called Taxonomy4_label
otu.table.005.nodes[,dim(otu.table.005.nodes)[2]+1] <- NA
colnames(otu.table.005.nodes)[dim(otu.table.005.nodes)[2]] <- "Taxonomy4_label"
# Fill in using loop- 
for(j in 1:length(x)) {
  for(i in 1:dim(otu.table.005.nodes)[1]) {
    if(is.na(otu.table.005.nodes$Taxonomy4)[i]) next
    if(otu.table.005.nodes$Taxonomy4[i] == x_labs[j]){otu.table.005.nodes$Taxonomy4_label[i] <- x[j]}
  }
}

# Next Taxonomy_5
rm(x,i,j)
x_labs <- as.array(unique(otu.table.005.nodes$Taxonomy5[!is.na(otu.table.005.nodes$Taxonomy5)])) # 
x <- 1:as.numeric(count(data.frame(x_labs))) # there are 187 categories
# Make new column for number index of Taxonomy5 called Taxonomy5_label
otu.table.005.nodes[,dim(otu.table.005.nodes)[2]+1] <- NA
colnames(otu.table.005.nodes)[dim(otu.table.005.nodes)[2]] <- "Taxonomy5_label"
# Fill in using loop- 
for(j in 1:length(x)) {
  for(i in 1:dim(otu.table.005.nodes)[1]) {
    if(is.na(otu.table.005.nodes$Taxonomy5)[i]) next
    if(otu.table.005.nodes$Taxonomy5[i] == x_labs[j]){otu.table.005.nodes$Taxonomy5_label[i] <- x[j]}
  }
}

# Next Taxonomy_6 
rm(x,i,j)
x_labs <- as.array(unique(otu.table.005.nodes$Taxonomy6[!is.na(otu.table.005.nodes$Taxonomy6)])) # 
x <- 1:as.numeric(count(data.frame(x_labs))) # there are 241 categories
# Make new column for number index of Taxonomy6 called Taxonomy6_label
otu.table.005.nodes[,dim(otu.table.005.nodes)[2]+1] <- NA
colnames(otu.table.005.nodes)[dim(otu.table.005.nodes)[2]] <- "Taxonomy6_label"
# Fill in using loop- 
for(j in 1:length(x)) {
  for(i in 1:dim(otu.table.005.nodes)[1]) {
    if(is.na(otu.table.005.nodes$Taxonomy6)[i]) next
    if(otu.table.005.nodes$Taxonomy6[i] == x_labs[j]){otu.table.005.nodes$Taxonomy6_label[i] <- x[j]}
  }
}



# Make network object
otu.table.005.net <- graph.data.frame(otu.table.005.links.sig, otu.table.005.nodes, directed=F)
# check attributes of network object
otu.table.005.net
E(otu.table.005.net) # shows 545 edges
V(otu.table.005.net) # there are 376 vertices
E(otu.table.005.net)$type
# V(otu.table.005.net)$Taxonomy # May want to come back later and make a column with full taxonomy
V(otu.table.005.net)$Taxonomy2
V(otu.table.005.net)$Taxonomy2_label

# Plot
# plot(otu.table.005.net) 

# Personalize plot
# Generate colors based Taxonomy2. Use the brewer palletes
unique(V(otu.table.005.net)$Taxonomy1) # there are 4 unique domains (including "No blast hit"), and also "NA"
# display.brewer.all()
# grDevices::colors
colrs <- c(brewer.pal(5, "Dark2"))
V(otu.table.005.net)$color <- colrs[V(otu.table.005.net)$Taxonomy1_label]

# Compute node degrees (#links) and use that to set node size: 
deg <- igraph::degree(otu.table.005.net, mode="all") # had to specifiy igraph::degree because there is another function degree
V(otu.table.005.net)$size <- 2*deg^0.9 # play with multiplication factors

# Turn off labels for now
V(otu.table.005.net)$label <- NA

# Set edge width based on weight: 
E(otu.table.005.net)$width <- E(otu.table.005.net)$weight*10

#change edge color: 
E(otu.table.005.net)$edge.color <- "black" 
#E(euk.net)$width <- 1+E(euk.net)$weight/12 

# plot
plot(otu.table.005.net) 
legend(x=-1.5, y=-1.1, c(unique(V(otu.table.005.net)$Taxonomy1)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
# Messed up because of the NAs but I am going to skip for now


# # Write a loop to generate all layouts and pick best ones
layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE) 
# # Remove layouts that do not apply to our graph. 
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]
# par(mfrow=c(3,3))
# for (layout in layouts) { 
#   print(layout) 
#   l <- do.call(layout, list(otu.table.005.net)) 
#   plot(otu.table.005.net, edge.arrow.mode=0, layout=l, main=layout)
# }
# 
# dev.off()

## Rendering all those plots takes a long time. Don't do it every time script is run.

l <- layout.circle(otu.table.005.net) 
plot(otu.table.005.net, layout=l) #vertex.label=V(otu.table.005.net))
legend(x=-1.5, y=-1.1, c(unique(V(otu.table.005.net)$Taxonomy1)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

# From this layout, you can see OTU #49 has most connection and is prominent in network
V(otu.table.005.net)$[409] # denovo ID = denovo100708 -E
otu.table.005.nodes[409,] # This is a syndiniales group II

# Can simplify this a bit by removing those vertices with few edges
cut.off <- .02*1^1.8 # = the corresponding size of a vertex(node) with only 1 edge
otu.table.005.net.mod <- delete_vertices(otu.table.005.net, V(otu.table.005.net)$size<cut.off)

E(otu.table.005.net.mod) # still shows 13810 edges
V(otu.table.005.net.mod) # and still 1327 vertices

l <- layout.auto(otu.table.005.net.mod) 
plot(otu.table.005.net.mod, layout=l)

# Next try plotting only positive or negative relationships
otu.table.005.net.pos <- otu.table.005.net - E(otu.table.005.net)[E(otu.table.005.net)$type=="negative"] 
otu.table.005.net.neg <- otu.table.005.net - E(otu.table.005.net)[E(otu.table.005.net)$type=="positive"] 

# plot(otu.table.005.net.pos, layout=l, main = "Positive Correlations", vertex.label=V(otu.table.005.net.pos)$Taxonomy3)
l <- layout.auto(otu.table.005.net) 
plot(otu.table.005.net.pos, layout=l, main = "Positive Correlations", vertex.label.font=20, vertex.label.cex = 2)
legend(x=-1.5, y=-1.1, c(unique(V(otu.table.005.net)$Taxonomy1)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

plot(otu.table.005.net.neg, layout=l, main = "Negative Correlations", vertex.label=V(otu.table.005.net.neg))
legend(x=-1.5, y=-1.1, c(unique(V(otu.table.005.net)$Taxonomy1)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)



# Export graph object to Gephi. Easier to manipulate/ plot 
write_graph(otu.table.005.net, "Network_PR2_OTUs/otu.table.005.net.gml", format = "gml") 


### In Gephi, I can plot just the neighbors of selected nodes but I can't view the reingold.tilford layout
# Look at reingold.tilford here, determine the major nodes, and use Gephi to show what's around those major nodes
l <- layout.reingold.tilford(otu.table.005.net) 
plot(otu.table.005.net, layout=l, vertex.label=V(otu.table.005.net))
legend(x=-1.5, y=-1.1, c(unique(V(otu.table.005.net)$Taxonomy1)), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

# From this layout, you can see OTU #49 has most connection and is most prominent in network
V(otu.table.005.net)$[409] # denovo ID = denovo100708 -E
otu.table.005.nodes[409,] # This is a syndiniales group II

# On a second level, nodes 19, 183, 364, 385, 392, 406, 506, 560, 587, 592, 677, 688, 792, 921, 963, 966, 985, 106, ...






#####  Make individual Networks for each redox regime ####

# Based on analysis from Suter et al 2018, can divide samples by oxycline, shallow anoxic, euxinic
# Take OTUs andenv variables from each, make new tables, re-run sparCC, and export to Gephi

# Use otu.table which already has all 3 domains (all OTUs, not filtered yet) 
# Split by redox regime
# Filter by abundance within that group of samples (so we don't have low abundance oxic OTUs in euxinic dataset, for example)
# Add metadeata
# abundance <0.5%, and cleaned so that any sample that doesn't have correpsonding data in other domains was removed

# Oxycline
otu.table.oxy <- otu.table[,c("OTU_ID", "A3a103A", "A3b103A",	"A3b198A",	"A3a234A","A2a200A", "A2a200B",
                                      "A3b234A", "A3b103B",	"A3a198B",	"A3b198B",	"A3a234B",	"A3b234B", 
                                      "A2a143A", "A2b143A",	"A2b200A",	"A2a237A",	"A2b237A", "A3a198A",
                                      "A2a143B", "A2b143B",	"A2b200B",	"A2a237B",	"A2b237B", "Taxonomy1",
                                      "Taxonomy2","Taxonomy3","Taxonomy4","Taxonomy5","Taxonomy6","Taxonomy7",
                                      "Taxonomy8","Taxonomy9")]

# Shallow Anoxic , remove A2b267A
otu.table.shallanox <- otu.table[,c("OTU_ID", "A3a295A",	"A3b295A",	"A3a314A", "A3b314A", "A2a267A",
                                            "A3a295B",	"A3b295B",	"A3a314B",	"A3b314B", "A2a247A",
                                            "A2b247A",	"A2b247B", "A2b267B", "A2a247B", "A2a267B", "Taxonomy1",
                                            "Taxonomy2","Taxonomy3","Taxonomy4","Taxonomy5","Taxonomy6","Taxonomy7",
                                            "Taxonomy8","Taxonomy9")]

# Euxinic, 
otu.table.eux <- otu.table[,c("OTU_ID", "A1b900A", "A3a900B",	"A1b900B", "A2a900AN", "A3a900A","A2b900AN", 
                              "A2a900BN", "A2b900BN","Taxonomy1","Taxonomy2","Taxonomy3","Taxonomy4","Taxonomy5",
                              "Taxonomy6","Taxonomy7","Taxonomy8","Taxonomy9")]

# Filter by most abundant. Retain any OTU that is at least 0.1% (RA = 0.001) in any sample
# Oxycline, reduced to 1233 OTUs
otu.table.oxy.001 <- filter(otu.table.oxy, A3a103A >= 0.001 | A3b103A >= 0.001 | A3b198A >= 0.001 | A3a234A >= 0.001 | 
                              A2a200A >= 0.001 | A2a200B >= 0.001 | A3b234A >= 0.001 | A3b103B >= 0.001 | A3a198B >= 0.001 | 
                              A3b198B >= 0.001 | A3a234B >= 0.001 | A3b234B >= 0.001 | A2a143A >= 0.001 | A2b143A >= 0.001 | 
                              A2b200A >= 0.001 | A2a237A >= 0.001 | A2b237A >= 0.001 | A3a198A >= 0.001 | A2a143B >= 0.001 | 
                              A2b143B >= 0.001 | A2b200B >= 0.001 | A2a237B >= 0.001 | A2b237B >= 0.001)

# Shallow anoxic, reduced to 586 OTUs
otu.table.shallanox.001 <- filter(otu.table.shallanox, A3a295A >= 0.001 | A3b295A >= 0.001 | A3a314A >= 0.001 | 
                                    A3b314A >= 0.001 | A2a267A >= 0.001 | A3a295B >= 0.001 | A3b295B >= 0.001 | 
                                    A3a314B >= 0.001 | A3b314B >= 0.001 | A2a247A >= 0.001 | A2b247A >= 0.001 | 
                                    A2b247B >= 0.001 | A2b267B >= 0.001 | A2a247B >= 0.001 | A2a267B >= 0.001)

# Euxinic, reduced to 547 OTUs
otu.table.eux.001 <- filter(otu.table.eux, A1b900A >= 0.001 | A3a900B >= 0.001 | A1b900B >= 0.001 | A2a900AN >= 0.001 |
                              A3a900A >= 0.001 | A2b900AN >= 0.001 | A2a900BN >= 0.001 | A2b900BN >= 0.001)



# Add in environmental data as rows
# upload metadata
meta <- read_delim("Metadata.txt",delim = "\t", col_names = T)

# Need to extract numeric metadata only
meta2 <- meta[6:23,]
meta2[,2:49] <- lapply(meta2[,2:49], as.numeric)

# Rename "key" to "OTU_ID" even though it's not really OTU- this will help match the dataframes
colnames(meta2)[1] <- colnames(otu.table)[1]

# Match metadata to OTU tables
otu.table.oxy.001 <- full_join(otu.table.oxy.001, meta2)
otu.table.oxy.001 <- otu.table.oxy.001[-c(34:58)]

otu.table.shallanox.001 <- full_join(otu.table.shallanox.001, meta2)
otu.table.shallanox.001 <- otu.table.shallanox.001[-c(6)]
otu.table.shallanox.001 <- otu.table.shallanox.001[-c(14)]
otu.table.shallanox.001 <- otu.table.shallanox.001[-c(24:56)]

otu.table.eux.001 <- full_join(otu.table.eux.001, meta2)
otu.table.eux.001 <- otu.table.eux.001[-c(7:9)]
otu.table.eux.001 <- otu.table.eux.001[-c(16:55)]


#  More filtering: remove samples for which there wasn't data for all 3 domains. Follow "SampleKey" excel sheet
# Oxycline: A3a198A, A2a200A, A2a200B
# Shall Anox: A3b314A, A2a267A, A2b267A, A2a247B, A2a267B
# Euxinic: A3a900A, A2b900AN, A2a900BN, A2b900BN
clean <- which(colnames(otu.table.oxy.001) %in% c("A3a198A", "A2a200A", "A2a200B"))
otu.table.oxy.001 <- otu.table.oxy.001[, -clean]

clean <- which(colnames(otu.table.shallanox.001) %in% c("A3b314A", "A2a267A", "A2b267A", "A2a247B", "A2a267B"))
otu.table.shallanox.001 <- otu.table.shallanox.001[, -clean]

clean <- which(colnames(otu.table.eux.001) %in% c("A3a900A", "A2b900AN", "A2a900BN", "A2b900BN"))
otu.table.eux.001 <- otu.table.eux.001[, -clean]

# Take note of which domains are which rows for each matrix
# Oxy matrix:
  # Rows 1-119 archaea
  # Rows 120-335 bacteria
  # Rows 336-1233 eukarya
  # 1234-1251 env var
# Shall Anox matrix:
  # Rows 1-61 archaea
  # Rows 62-183 bacteria (although last 2 OTUs have no blast hit at domain level, they are from bac dataset)
  # Rows 184-586 eukarya
  # Rows 587-604 env var
# Euxinic matrix
  # Rows 1-85 archaea (one OTU is also "No blast hit"at domain level")
  # Rows 86-190 (2 OTUs "Not blast hit" at domain level)
  # Rows 191-547 eukarya
  # Rows 548-565 env var


# Export for SparCC
write.table(otu.table.oxy.001[,1:21], "Network_PR2_OTUs/otu.table.oxy.001.txt", sep = "\t",row.names = FALSE, quote=FALSE)
write.table(otu.table.shallanox.001[,12], "Network_PR2_OTUs/otu.table.shallanox.001.txt", sep = "\t",row.names = FALSE, quote=FALSE)
write.table(otu.table.eux.001[,1:5], "Network_PR2_OTUs/otu.table.eux.001.txt", sep = "\t",row.names = FALSE, quote=FALSE)

# Run correlation analyses on these in SparCC in Python
# See lab notebook for code







# Import results from SparCC and make network objects
# Oxycline
otu.table.oxy.001.corr <- as.matrix(read.table("Network_PR2_OTUs/otu.table.oxy.001_basis_corr/cor_sparcc.out", header = T, row.names = 1, sep = "\t"))
otu.table.oxy.001.oneside <- as.matrix(read.table("Network_PR2_OTUs/otu.table.oxy.001_pvals/pvals_one_sided.txt", header = T, row.names = 1, sep = "\t"))

# Shallow Anoxic
otu.table.shallanox.001.corr <- as.matrix(read.table("Network_PR2_OTUs/otu.table.shallanox.001_basis_corr/cor_sparcc.out", header = T, row.names = 1, sep = "\t"))
otu.table.shallanox.001.oneside <- as.matrix(read.table("Network_PR2_OTUs/otu.table.shallanox.001_pvals/pvals_one_sided.txt", header = T, row.names = 1, sep = "\t"))

# Euxinic
otu.table.eux.001.corr <- as.matrix(read.table("Network_PR2_OTUs/otu.table.eux.001_basis_corr/cor_sparcc.out", header = T, row.names = 1, sep = "\t"))
otu.table.eux.001.oneside <- as.matrix(read.table("Network_PR2_OTUs/otu.table.eux.001_pvals/pvals_one_sided.txt", header = T, row.names = 1, sep = "\t"))





# Remove everything above the diaganol in each matrix so results are not repeated later on when these are linearized
otu.table.oxy.001.corr[upper.tri(otu.table.oxy.001.corr)] <- NA
otu.table.oxy.001.oneside[upper.tri(otu.table.oxy.001.oneside)] <- NA

otu.table.shallanox.001.corr[upper.tri(otu.table.shallanox.001.corr)] <- NA
otu.table.shallanox.001.oneside[upper.tri(otu.table.shallanox.001.oneside)] <- NA

otu.table.eux.001.corr[upper.tri(otu.table.eux.001.corr)] <- NA
otu.table.eux.001.oneside[upper.tri(otu.table.eux.001.oneside)] <- NA




# Associate the denovo IDs with a marker that indicates which of the 3 domains (because there are denovo repeats when datasets are combined)
# From above-
# Oxy matrix:
    # Rows 1-119 archaea
    # Rows 120-335 bacteria
    # Rows 336-1233 eukarya
    # 1234-1251 env var
colnames(otu.table.oxy.001.corr)[1:119] <- paste(colnames(otu.table.oxy.001.corr)[1:119],"-A")
rownames(otu.table.oxy.001.corr)[1:119] <- paste(rownames(otu.table.oxy.001.corr)[1:119],"-A")
colnames(otu.table.oxy.001.corr)[120:335] <- paste(colnames(otu.table.oxy.001.corr)[120:335],"-B")
rownames(otu.table.oxy.001.corr)[120:335] <- paste(rownames(otu.table.oxy.001.corr)[120:335],"-B")
colnames(otu.table.oxy.001.corr)[336:1233] <- paste(colnames(otu.table.oxy.001.corr)[336:1233],"-E")
rownames(otu.table.oxy.001.corr)[336:1233] <- paste(rownames(otu.table.oxy.001.corr)[336:1233],"-E")

colnames(otu.table.oxy.001.oneside)[1:119] <- paste(colnames(otu.table.oxy.001.oneside)[1:119],"-A")
rownames(otu.table.oxy.001.oneside)[1:119] <- paste(rownames(otu.table.oxy.001.oneside)[1:119],"-A")
colnames(otu.table.oxy.001.oneside)[120:335] <- paste(colnames(otu.table.oxy.001.oneside)[120:335],"-B")
rownames(otu.table.oxy.001.oneside)[120:335] <- paste(rownames(otu.table.oxy.001.oneside)[120:335],"-B")
colnames(otu.table.oxy.001.oneside)[336:1233] <- paste(colnames(otu.table.oxy.001.oneside)[336:1233],"-E")
rownames(otu.table.oxy.001.oneside)[336:1233] <- paste(rownames(otu.table.oxy.001.oneside)[336:1233],"-E")

# Shall Anox matrix:
    # Rows 1-61 archaea
    # Rows 62-183 bacteria (although last 2 OTUs have no blast hit at domain level, they are from bac dataset)
    # Rows 184-586 eukarya
    # Rows 587-604 env var
colnames(otu.table.shallanox.001.corr)[1:61] <- paste(colnames(otu.table.shallanox.001.corr)[1:61],"-A")
rownames(otu.table.shallanox.001.corr)[1:61] <- paste(rownames(otu.table.shallanox.001.corr)[1:61],"-A")
colnames(otu.table.shallanox.001.corr)[62:183] <- paste(colnames(otu.table.shallanox.001.corr)[62:183],"-B")
rownames(otu.table.shallanox.001.corr)[62:183] <- paste(rownames(otu.table.shallanox.001.corr)[62:183],"-B")
colnames(otu.table.shallanox.001.corr)[184:586] <- paste(colnames(otu.table.shallanox.001.corr)[184:586],"-E")
rownames(otu.table.shallanox.001.corr)[184:586] <- paste(rownames(otu.table.shallanox.001.corr)[184:586],"-E")

colnames(otu.table.shallanox.001.oneside)[1:61] <- paste(colnames(otu.table.shallanox.001.oneside)[1:61],"-A")
rownames(otu.table.shallanox.001.oneside)[1:61] <- paste(rownames(otu.table.shallanox.001.oneside)[1:61],"-A")
colnames(otu.table.shallanox.001.oneside)[62:183] <- paste(colnames(otu.table.shallanox.001.oneside)[62:183],"-B")
rownames(otu.table.shallanox.001.oneside)[62:183] <- paste(rownames(otu.table.shallanox.001.oneside)[62:183],"-B")
colnames(otu.table.shallanox.001.oneside)[184:586] <- paste(colnames(otu.table.shallanox.001.oneside)[184:586],"-E")
rownames(otu.table.shallanox.001.oneside)[184:586] <- paste(rownames(otu.table.shallanox.001.oneside)[184:586],"-E")

# Euxinic matrix
    # Rows 1-85 archaea (one OTU is also "No blast hit"at domain level")
    # Rows 86-190 (2 OTUs "Not blast hit" at domain level)
    # Rows 191-547 eukarya
    # Rows 548-565 env var
colnames(otu.table.eux.001.corr)[1:85] <- paste(colnames(otu.table.eux.001.corr)[1:85],"-A")
rownames(otu.table.eux.001.corr)[1:85] <- paste(rownames(otu.table.eux.001.corr)[1:85],"-A")
colnames(otu.table.eux.001.corr)[86:190] <- paste(colnames(otu.table.eux.001.corr)[86:190],"-B")
rownames(otu.table.eux.001.corr)[86:190] <- paste(rownames(otu.table.eux.001.corr)[86:190],"-B")
colnames(otu.table.eux.001.corr)[191:547] <- paste(colnames(otu.table.eux.001.corr)[191:547],"-E")
rownames(otu.table.eux.001.corr)[191:547] <- paste(rownames(otu.table.eux.001.corr)[191:547],"-E")

colnames(otu.table.eux.001.oneside)[1:85] <- paste(colnames(otu.table.eux.001.oneside)[1:85],"-A")
rownames(otu.table.eux.001.oneside)[1:85] <- paste(rownames(otu.table.eux.001.oneside)[1:85],"-A")
colnames(otu.table.eux.001.oneside)[86:190] <- paste(colnames(otu.table.eux.001.oneside)[86:190],"-B")
rownames(otu.table.eux.001.oneside)[86:190] <- paste(rownames(otu.table.eux.001.oneside)[86:190],"-B")
colnames(otu.table.eux.001.oneside)[191:547] <- paste(colnames(otu.table.eux.001.oneside)[191:547],"-E")
rownames(otu.table.eux.001.oneside)[191:547] <- paste(rownames(otu.table.eux.001.oneside)[191:547],"-E")






# Make a "links" matrix (the edges) with correlation weights of all pairs of nodes
# Simplify correlation matrix into a list of pairs
otu.table.oxy.001.links <- melt(otu.table.oxy.001.corr)
colnames(otu.table.oxy.001.links)[3] <- "weight"
colnames(otu.table.oxy.001.links)

otu.table.shallanox.001.links <- melt(otu.table.shallanox.001.corr)
colnames(otu.table.shallanox.001.links)[3] <- "weight"
colnames(otu.table.shallanox.001.links)

otu.table.eux.001.links <- melt(otu.table.eux.001.corr)
colnames(otu.table.eux.001.links)[3] <- "weight"
colnames(otu.table.eux.001.links)




# Put pvalues in this list for now
otu.table.oxy.001.pvals <- melt(otu.table.oxy.001.oneside)
otu.table.oxy.001.links <- inner_join(otu.table.oxy.001.links, otu.table.oxy.001.pvals)
colnames(otu.table.oxy.001.links)[1] <- "from"
colnames(otu.table.oxy.001.links)[2] <- "to"
colnames(otu.table.oxy.001.links)[4] <- "pval_1side"
colnames(otu.table.oxy.001.links)

otu.table.shallanox.001.pvals <- melt(otu.table.shallanox.001.oneside)
otu.table.shallanox.001.links <- inner_join(otu.table.shallanox.001.links, otu.table.shallanox.001.pvals)
colnames(otu.table.shallanox.001.links)[1] <- "from"
colnames(otu.table.shallanox.001.links)[2] <- "to"
colnames(otu.table.shallanox.001.links)[4] <- "pval_1side"
colnames(otu.table.shallanox.001.links)

otu.table.eux.001.pvals <- melt(otu.table.eux.001.oneside)
otu.table.eux.001.links <- inner_join(otu.table.eux.001.links, otu.table.eux.001.pvals)
colnames(otu.table.eux.001.links)[1] <- "from"
colnames(otu.table.eux.001.links)[2] <- "to"
colnames(otu.table.eux.001.links)[4] <- "pval_1side"
colnames(otu.table.eux.001.links)


# Filter out non signficant relationships. then remove "pval" column
# Use p-value cutoff of 0.05 
pval <- 0.05

otu.table.oxy.001.links.sig <- filter(otu.table.oxy.001.links, pval_1side <=pval)
otu.table.oxy.001.links.sig <- otu.table.oxy.001.links.sig[-c(4)]
# Reduced from 1,565,001 links to 93,184 

otu.table.shallanox.001.links.sig <- filter(otu.table.shallanox.001.links, pval_1side <=pval)
otu.table.shallanox.001.links.sig <- otu.table.shallanox.001.links.sig[-c(4)]
# Reduced from 364,186 links to 21,679

otu.table.eux.001.links.sig <- filter(otu.table.eux.001.links, pval_1side <=pval)
otu.table.eux.001.links.sig <- otu.table.eux.001.links.sig[-c(4)]
# Reduced from 319,225 links to 19,175




# Add in "type" column to show positive or negative relationship. Then just retain the absolute value of the correlation
otu.table.oxy.001.links.sig$type <- ifelse(otu.table.oxy.001.links.sig$weight > 0, 'positive',
                                                 ifelse(otu.table.oxy.001.links.sig$weight < 0, 'negative',"NA")
)
otu.table.oxy.001.links.sig$weight <- abs(otu.table.oxy.001.links.sig$weight)


otu.table.shallanox.001.links.sig$type <- ifelse(otu.table.shallanox.001.links.sig$weight > 0, 'positive',
                                       ifelse(otu.table.shallanox.001.links.sig$weight < 0, 'negative',"NA")
)
otu.table.shallanox.001.links.sig$weight <- abs(otu.table.shallanox.001.links.sig$weight)


otu.table.eux.001.links.sig$type <- ifelse(otu.table.eux.001.links.sig$weight > 0, 'positive',
                                           ifelse(otu.table.eux.001.links.sig$weight < 0, 'negative',"NA")
)
otu.table.eux.001.links.sig$weight <- abs(otu.table.eux.001.links.sig$weight)







# Remove those with correl coef less than 0.25
coef <- 0.2

otu.table.oxy.001.links.sig <- filter(otu.table.oxy.001.links.sig, weight >= coef)
# Reduced from 93184 links to 1477

otu.table.shallanox.001.links.sig <- filter(otu.table.shallanox.001.links.sig, weight >= coef)
# Reduced from 21679 links to 5483

otu.table.eux.001.links.sig <- filter(otu.table.eux.001.links.sig, weight >= coef)
# Reduced from 19175 links to 19174... so everything here is highly correlated but remember sample size is smaller





# # Next remove the edges that are connected env variables- these are not interesting for my purposes and they force the plot a certain way that is not meaningful
# otu.table.shallanox.001.links.sig <- otu.table.shallanox.001.links.sig[-c(1148:1167),]


# Redefine the vertices (nodes) by just what's left in the matrix after trimming out non sig relationships
x_oxy <- melt(t(otu.table.oxy.001.links.sig[1:2]))
x_oxy <- data.frame(unique(x_oxy$value))
count(x_oxy) # there are 1146 unique OTUs still remaining
colnames(x_oxy) <- "OTU_ID"

x_shallanox <- melt(t(otu.table.shallanox.001.links.sig[1:2]))
x_shallanox <- data.frame(unique(x_shallanox$value))
count(x_shallanox) # there are 604 unique OTUs still remaining
colnames(x_shallanox) <- "OTU_ID"

x_eux <- melt(t(otu.table.eux.001.links.sig[1:2]))
x_eux <- data.frame(unique(x_eux$value))
count(x_eux) # there are 565 unique OTUs still remaining
colnames(x_eux) <- "OTU_ID"




# Join the retained OTU IDs to taxonomy- this is the new "nodes" dataset

# I think the easiest way to do this is take the taxonomy/ OTU ID from arch_abun_table, add in "-A" 
# to the denovo IDs, then match by left_join
arch_taxonomy <- arch.otu.table[-c(2:43)]
arch_taxonomy$OTU_ID <- paste(arch_taxonomy$OTU_ID,"-A")
bac_taxonomy <- bac.otu.table[-c(2:46)]
bac_taxonomy$OTU_ID <- paste(bac_taxonomy$OTU_ID,"-B")
euk_taxonomy <- euk.otu.table[-c(2:48)]
euk_taxonomy$OTU_ID <- paste(euk_taxonomy$OTU_ID,"-E")

taxonomy <- full_join(arch_taxonomy, bac_taxonomy) #721776 rows, 11 cols
taxonomy <- full_join(taxonomy, euk_taxonomy) # 1167818, 11 cols

# Then use this to make nodes for each network
otu.table.oxy.001.nodes <- left_join(x_oxy, taxonomy)
otu.table.shallanox.001.nodes <- left_join(x_shallanox, taxonomy)
otu.table.eux.001.nodes <- left_join(x_eux, taxonomy)





# Make network objects
otu.table.oxy.001.net <- graph.data.frame(otu.table.oxy.001.links.sig, otu.table.oxy.001.nodes, directed=F)
otu.table.shallanox.001.net <- graph.data.frame(otu.table.shallanox.001.links.sig, otu.table.shallanox.001.nodes, directed=F)
otu.table.eux.001.net <- graph.data.frame(otu.table.eux.001.links.sig, otu.table.eux.001.nodes, directed=F)

# Export graph object to Gephi. Easier to manipulate/ plot 
write_graph(otu.table.oxy.001.net, "Network_PR2_OTUs/otu.table.oxy.001.net.gml", format = "gml") 
write_graph(otu.table.shallanox.001.net, "Network_PR2_OTUs/otu.table.shallanox.001.net.gml", format = "gml") 
write_graph(otu.table.eux.001.net, "Network_PR2_OTUs/otu.table.eux.001.net.gml", format = "gml") 











# Calculate some network properties using igraph :

# Mean Clustering Coefficient (Transitivity) of Network
transitivity(otu.table.oxy.001.net, type = c("undirected"))
transitivity(otu.table.shallanox.001.net, type = c("undirected"))
transitivity(otu.table.eux.001.net, type = c("undirected"))




# Degree
mean(igraph::degree(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), mode = c("all")))
mean(igraph::degree(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), mode = c("all")))
mean(igraph::degree(otu.table.eux.001.net, v = V(otu.table.eux.001.net), mode = c("all")))

# Pick out nodes with highest degree
# Oxy
sort(igraph::degree(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), mode = c("all")))
# Degrees range from 1 to 4. Pick out all 3s and 4s
oxy_hi_deg <- which(igraph::degree(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), mode = c("all")) >3)

# ShallAnox
sort(igraph::degree(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), mode = c("all")))
# Degrees range from 2 to 34. Take out highest values, anything greater than 20
shallanox_hi_deg <- which(igraph::degree(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), mode = c("all")) > 20)

# Euxinic
sort(igraph::degree(otu.table.eux.001.net, v = V(otu.table.eux.001.net), mode = c("all")))
# Degrees range from 39 to 96. Take out highest values, anything greater than 70
eux_hi_deg <- which(igraph::degree(otu.table.eux.001.net, v = V(otu.table.eux.001.net), mode = c("all")) > 70)




# Avg Shortest Path of Network
mean_distance(otu.table.oxy.001.net) 
mean_distance(otu.table.shallanox.001.net) 
mean_distance(otu.table.eux.001.net) 




# Avg Betweenness centrality of Network
mean(igraph::betweenness(otu.table.oxy.001.net))
mean(igraph::betweenness(otu.table.shallanox.001.net))
mean(igraph::betweenness(otu.table.eux.001.net))

# Pick out nodes with lowest betweenness values--> these are candidates for keystone species
# Oxy
y <- which(igraph::betweenness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), directed = FALSE) == min(igraph::betweenness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), directed = FALSE)))
igraph::betweenness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), directed = FALSE)[y]
# List of OTUs with low values
oxy_low_btwnss <- sort(igraph::betweenness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net), directed = FALSE))[1:294]
# I considered the bottom 294 for this network since the 294th has betweenness value of 2 and the 295th OTU 
#  has  value of 335

# Shallanox
# List of OTUs with the minimum value (only 1 OTU with value of 52 in this case)
z <- which(igraph::betweenness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), directed = FALSE, weights =) == min(igraph::betweenness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), directed = FALSE)))
igraph::betweenness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), directed = FALSE)[z]
# List of OTUs with low values
shallanox_low_btwnss <- sort(igraph::betweenness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net), directed = FALSE))[1:294]
# For this network, the betweenness values of OTUs climbs very steadily, so just take bottom 294

# Euxinic
# List of OTUs with the minimum value (only 2 OTUs with value of 72 in this case)
w <- which(igraph::betweenness(otu.table.eux.001.net, v = V(otu.table.eux.001.net), directed = FALSE) == min(igraph::betweenness(otu.table.eux.001.net, v = V(otu.table.eux.001.net), directed = FALSE)))
igraph::betweenness(otu.table.eux.001.net, v = V(otu.table.eux.001.net), directed = FALSE)[w]
# List of OTUs with low values
eux_low_btwnss <- sort(igraph::betweenness(otu.table.eux.001.net, v = V(otu.table.eux.001.net), directed = FALSE))[1:294]
# For this network, the betweenness values of OTUs climbs very steadily, so just take bottom 294




# Avg Closeness centrality of Network
mean(igraph::closeness(otu.table.oxy.001.net)) # get warning about this one:"closeness centrality is not well-defined for disconnected graphs"
mean(igraph::closeness(otu.table.shallanox.001.net))
mean(igraph::closeness(otu.table.eux.001.net))

# Pick out nodes with highest closeness values--> these are candidates for keystone species
# Oxy
sort(igraph::closeness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net)))
max(igraph::closeness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net)))
oxy_hi_closnss <- which(igraph::closeness(otu.table.oxy.001.net, v = V(otu.table.oxy.001.net)) > 2.21e-05)

# ShallAnox
sort(igraph::closeness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net)))
max(igraph::closeness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net)))
shallanox_hi_closnss <- which(igraph::closeness(otu.table.shallanox.001.net, v = V(otu.table.shallanox.001.net)) > 0.003)
# Note that one of these, denovo207699, is the No Blast hit OTU 
# I blasted the sequence and the top hit is Thioalkalivibrio sulfidiphilus, an obligate S-oxidizing chemoautotroph
# that had its genome sequenced in 2017 (we were using SILVA annotations from before then). 
# But it's only an 84% match.

# Euxinic
sort(igraph::closeness(otu.table.eux.001.net, v = V(otu.table.eux.001.net)))
max(igraph::closeness(otu.table.eux.001.net, v = V(otu.table.eux.001.net)))
eux_hi_closnss <- which(igraph::closeness(otu.table.eux.001.net, v = V(otu.table.eux.001.net)) > 0.00265)



# Next, look for overlapp in OTUs between among the lists of low betweeness, low closeness, and high degree
# Oxy network
length(oxy_hi_deg)
length(oxy_low_btwnss) #longest list
length(oxy_hi_closnss) #shortest

match(names(oxy_hi_deg), names(oxy_hi_closnss)) 
# Translation: in position 6 and 7 of oxy_hi_deg there is a match to pos 1 and 2 from oxy_hi_closnss
# "denovo253297 -E"       "denovo447012 -E" 
# "denovo253297 -E"       "denovo447012 -E" 
# Extract those names and check for them in the third matrix, oxy_low_btwnss
oxy_temp <- names(oxy_hi_deg)[!is.na(match(names(oxy_hi_deg), names(oxy_hi_closnss)))]

match(names(oxy_low_btwnss), oxy_temp) 
# no matches, so no keystone candidates!



# Shallow Anoxic network
length(shallanox_hi_deg)
length(shallanox_low_btwnss) 
length(shallanox_hi_closnss) 

match(names(shallanox_hi_deg), names(shallanox_low_btwnss)) 
shallanox_temp <- names(shallanox_hi_deg)[!is.na(match(names(shallanox_hi_deg), names(shallanox_low_btwnss)))]

match(names(shallanox_hi_closnss), shallanox_temp) 
shallanox_keystone_cands <- names(shallanox_hi_closnss)[!is.na(match(names(shallanox_hi_closnss), shallanox_temp))]
# There are 7 keystone candidates. One is NO3, so there are 6 keystone candidate OTUs:
# "denovo158333 -B" "denovo64710 -E"  "denovo431262 -E" "denovo244612 -B" "denovo49812 -E"  "denovo294767 -B"




# Shallow Anoxic network
length(eux_hi_deg)
length(eux_low_btwnss) 
length(eux_hi_closnss) 

match(names(eux_hi_closnss), names(eux_hi_deg)) 
eux_temp <- names(eux_hi_closnss)[!is.na(match(names(eux_hi_closnss), names(eux_hi_deg)))]

match(names(eux_low_btwnss), eux_temp) 
eux_keystone_cands <- names(eux_low_btwnss)[!is.na(match(names(eux_low_btwnss), eux_temp))]
# There are 7 keystone candidates:
# "denovo21694 -E"  "denovo92288 -A"  "denovo48360 -A"  "denovo161618 -E" "denovo288601 -E" "denovo106090 -E" "denovo165566 -A"



# Next: Make network figures highlightingthe dynamics of these keystones
# Also plot some abundances of Cariacotrichea, for example, to make guesses about its symbionts



##### 
test <- which(otu.table.shallanox$OTU_ID %in% c("denovo231149","denovo313884"))
# plot row 600530 (BS-GSO-2) vs row 984927 (Cariacotrichea)
BSGSO2 <- otu.table[600530,2:43]
Cariacotrichea <- otu.table[984927,2:43]
x <- c(600530, 984927)
temp <- otu.table[x,]
#plot(temp[1,2:43],temp[2,2:43])

temp[subset(temp,!is.na(temp[1,]))]

install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")

chart.Correlation(temp, histogram=TRUE, pch=19)

plot(subset$Time[!is.na(subset$A)],subset$A[!is.na(subset$A)],type="l")  



