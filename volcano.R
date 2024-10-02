# Set the folder containing current R script as working directory
setwd(here::here())
print(paste("Current working dir:", getwd()))

# Load R packages
lapply(c("EnhancedVolcano", "tidyverse", "org.Hs.eg.db"), 
       require, character.only = TRUE)


################################################################################

# Read in male RNA-seq data
data_male <- read.csv("./Differential_expression_analysis_table_male.csv")
uniq_name <- make.names(data_male$Gene.name, unique = TRUE)
row.names(data_male) <- uniq_name

values_to_exclude <- c('Rpl3-ps1', 'Dennd4b', 'Kif13b', '2-Mar', 'Kif13b', 'Pcca', 'A230057D06Rik', 'Crybg3', 'Deaf1', 'Dlg2', 'Gm15853', 'Gm16499', 'Gm16701', 'Gm2464', 'Gm4430', 'Ksr2', 'Magi1', 'Nudt8', 'Pcdha11', 'Ppat', 'Rit1', 'Spag6', 'St6galnac2', 'U2af1l4', 'Zkscan7')
data_male <- data_male[!data_male$Gene.name %in% values_to_exclude, ]

# Read in female RNA-seq data
data_female <- read.csv("./Differential_expression_analysis_table_female.csv")
uniq_name <- make.names(data_female$Gene.name, unique = TRUE)
row.names(data_female) <- uniq_name

values_to_exclude <- c('mt-Tq', 'Gm37928')
data_female <- data_female[!data_female$Gene.name %in% values_to_exclude, ]


################################################################################

# Plot male
p <- EnhancedVolcano(data_male,
                     title = 'ATF4 male (8 weeks) RNA-seq',
                     lab = rownames(data_male),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     xlim = c(-3,3),
                     ylim = c(0,55), 
                     pCutoff = 10e-5,
                     FCcutoff = 0.4,)
ggsave("ATF4_male_volcano.pdf", plot = p, width = 7, height = 7, units = "in", dpi = 1200)
rm("p")

# Plot female
p <- EnhancedVolcano(data_female,
                     title = 'ATF4 female (8 weeks) RNA-seq',
                     lab = rownames(data_female),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     xlim = c(-3,3),
                     ylim = c(0,32), 
                     pCutoff = 10e-5,
                     FCcutoff = 0.4,)
ggsave("ATF4_female_volcano.pdf", plot = p, width = 7, height = 7, units = "in", dpi = 1200)
rm("p")


################################################################################

# exit
rm(list = ls())
