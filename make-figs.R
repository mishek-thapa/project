library(tidyverse)
library(vtree)
library(flextable)
library(magick)
library(ggpubr)


# Code is organized by the section in which each figures and tables came from.
# The code for the eda takes a long time to install and run. 


###########
# EDA Figures
###########

#code from: https://stackoverflow.com/questions/57175351/flextable-autofit-in-a-rmarkdown-to-word-doc-causes-table-to-go-outside-page-mar

FitFlextableToPage <- function(ft, pgwidth = 6){
  
  ft_out <- ft %>% autofit()
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}

#idep functions: http://bioinformatics.sdstate.edu/idep/
source('iDEP_core_functions.R') 


# Input files 
# Expression file has to use Ensembl for gene ID. Otherwise, use custom pathway database with matching IDs. 
inputFile <- 'eda_data.csv'  # Expression matrix
# Experiment design file
sampleInfoFile <- 'cng_design_file.csv'  
#Gene symbols, location etc. 
geneInfoFile <- 'Mouse__mmusculus_gene_ensembl_GeneInfo.csv' 
# pathway database in SQL; can be GMT format 
geneSetFile <- 'Mouse__mmusculus_gene_ensembl.db'   
STRING10_speciesFile <- 'https://raw.githubusercontent.com/iDEP-SDSU/idep/master/shinyapps/idep/STRING10_species.csv' 

# Parameters
input_missingValue <- 'geneMedian'	#Missing values imputation method
input_dataFileFormat <- 1	#1- read counts, 2 FKPM/RPKM or DNA microarray
input_minCounts <- 0.5	#Min counts
input_NminSamples <- 1	#Minimum number of samples 
input_countsLogStart <- 4	#Pseudo count for log CPM
input_CountsTransform <- 3	#Methods for data transformation of counts. 1-EdgeR's logCPM; 2-VST; 3-rlog 

#Read data files
readData.out <- readData(inputFile) 
readSampleInfo.out <- readSampleInfo(sampleInfoFile) 
input_selectOrg ="NEW" 
input_selectGO <- NULL 	#Gene set category 
input_noIDConversion = TRUE  
allGeneInfo.out <- geneInfo(geneInfoFile) 
converted.out = NULL 
convertedData.out <- convertedData()	 
nGenesFilter()  
convertedCounts.out <- convertedCounts()  # converted counts, just for compatibility 
readCountsBias()  # detecting bias in sequencing depth


input_selectFactors <- 'genotype'	#Factor coded by color
#PCA plots


library(ggpubr)
pca_fig <-  PCAplot()	+ theme_pubr() + labs(color = "Genotype", 
                                            title = "Groups well separated across PCA") + 
  scale_color_brewer(palette = "Dark2",
                     labels = c("CNG w/ 30 Day Treatment",
                                "CNG",
                                "Control")) +
  theme(legend.position = "bottom")
dev.off()


##########################
# 4. k-Means clustering 
##########################
input_nGenesKNN <- 2000	#Number of genes fro k-Means
input_nClusters <- 4	#Number of clusters 
maxGeneClustering = 12000
input_kmeansNormalization <- 'geneMean'	#Normalization
input_KmeansReRun <- 0	#Random seed 

distributionSD()  #Distribution of standard deviations 
KmeansNclusters()  #Number of clusters 
Kmeans.out = Kmeans()   #Running K-means 



png(filename = "heatmap.png", width = 800, height = 800, pointsize = 20)
KmeansHeatmap()   #Heatmap for k-Means
dev.off()


###########
# Results Figures
###########

library(EBSeqHMM)

counts <- read_csv("final_all_counts.csv") %>%
  select(1, 5:12, 16:19)

GeneExampleData <- as.matrix(counts[-1])
row.names(GeneExampleData) = counts$gene

CondVector <- rep(paste("t",1:3,sep=""),each=4)
print(CondVector)



Conditions <- factor(CondVector, levels=c("t1","t2","t3"))
str(Conditions)
levels(Conditions)


Sizes <- MedianNorm(GeneExampleData)

GeneNormData <- GetNormalizedMat(GeneExampleData, Sizes)


PlotExp(GeneNormData, Conditions, Name="ENSMUSG00000025917")

EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=10)

GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.1)

head(GeneDECalls)
str(GeneDECalls)

GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=.1,cutoff=.6, OnlyDynamic=TRUE)
str(GeneConfCalls$EachPathNames)


dd <- unlist(GeneConfCalls$EachPathNames["Down-Down"])
uu <- unlist(GeneConfCalls$EachPathNames["Up-Up"])
du <- unlist(GeneConfCalls$EachPathNames["Down-Up"])
ud <- unlist(GeneConfCalls$EachPathNames["Up-Down"])



#plotgenes function

norm_counts <- as.data.frame(GeneNormData)
norm_counts$gene <- counts$gene

colnames(norm_counts)[-13] = CondVector

library(janitor)



d1 <- norm_counts %>%
  clean_names() %>%
  pivot_longer(cols = 1:4,
               names_to = "t1",
               values_to = "t1_count") %>%
  select(9,11)

d2 <- norm_counts %>%
  clean_names() %>%
  pivot_longer(cols = 5:8,
               names_to = "t2",
               values_to = "t2_count") %>%
  select(11)


d3 <- norm_counts %>%
  clean_names() %>%
  pivot_longer(cols = 9:12,
               names_to = "t3",
               values_to = "t3_count")  %>%
  select(11)


norm_df <- cbind(d1, d2, d3)

sig_genes <- norm_df %>%
  group_by(gene) %>%
  summarise(P30 = median(t1_count),
            P90 = median(t2_count),
            P210 = median(t3_count)) %>%
  mutate(signi = case_when(
    gene %in% dd ~ "Down-Down",
    gene %in% uu ~ "upup",
    gene %in% ud ~ "updown",
    gene %in% du ~ "Down-Up",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(signi)) %>%
  pivot_longer(cols = 2:4,
               names_to = "timepoint",
               values_to = "n_expression")

counts <- sig_genes %>%
  dplyr::count(signi) %>%
  mutate(label = as.character(n/3)
  )

sig_genes$timepoint <- factor(sig_genes$timepoint, levels = c("P30", "P90", "P210"))


sig_genes <- sig_genes %>%
  mutate(mygene = case_when(
    gene %in% c("ENSMUSG00000051228",
                "ENSMUSG00000000617",
                "ENSMUSG00000030523",
                "ENSMUSG00000070337",
                "ENSMUSG00000021719",
                "ENSMUSG00000024186",
                "ENSMUSG00000025739",
                "ENSMUSG00000004630",
                "ENSMUSG00000031748",
                "ENSMUSG00000023439",
                'ENSMUSG00000056043',
                'ENSMUSG00000032192',
                "ENSMUSG00000015968") ~ "Post-synaptic",
    gene %in% c("ENSMUSG00000093865",
                "ENSMUSG00000048988",
                "ENSMUSG00000043460",
                "ENSMUSG00000039952",
                "ENSMUSG00000042961",
                "ENSMUSG00000031142",
                "ENSMUSG00000045103") ~ "lrit3",
    gene %in% c("ENSMUSG00000024842,
                ENSMUSG00000041460") ~ "CSNB",
    TRUE ~ NA_character_))

mygenes <- sig_genes %>%
  filter(!is.na(mygene))



(fig1 <- sig_genes %>%
    filter(is.na(mygene)) %>%
    ggplot(aes(x = timepoint,
               y = log(n_expression), group = gene)) +
    theme_pubr() +
    geom_line(alpha = 0.2, color = "grey") +
    geom_point(alpha = 0.2, color = "grey") +
    facet_wrap(vars(signi), ncol =1, drop = TRUE, scales = "free_y") +
    geom_line(data = mygenes, aes(x = timepoint,
                                  y = log(n_expression), group = mygene, color = mygene)) +
    geom_point(data = mygenes, aes(x = timepoint,
                                   y = log(n_expression), group = mygene, color = mygene)) +
    geom_text(data    = counts,
              mapping = aes(x = Inf, y = Inf, label = label),
              hjust   = 1,
              vjust   = 1,
              inherit.aes=FALSE
    ) +  
    geom_text(data    = counts,
              mapping = aes(x = Inf, y = Inf, label = label),
              hjust   = 1,
              vjust   = 1,
              inherit.aes=FALSE
    ) +
    labs(title = "Cngb1",
         y = "Log of Normalized Expression",
         x = "Time Point",
         color = "")) 


down <- sig_genes %>%
  filter(signi %in% c("downdown"),
         timepoint %in% c("P30", "P210")) %>%
  pivot_wider(id_cols = gene, names_from  = timepoint, values_from = n_expression) %>%
  mutate(log2_fold_change = log2(P30/P210)) %>%
  arrange(desc(log2_fold_change)) %>%
  slice(1:10) %>%
  mutate(change = "downregulation")

up <- sig_genes %>%
  filter(signi %in% c("upup"),
         timepoint %in% c("P30", "P210")) %>%
  pivot_wider(id_cols = gene, names_from  = timepoint, values_from = n_expression) %>%
  mutate(log2_fold_change = log2(P30/P210)) %>%
  arrange(log2_fold_change, desc(P210)) %>%
  slice(1:10) %>%
  mutate(change = "upregulation")

cng<- rbind(up, down)
cng$condition <- "cng"






# control analysis 

counts <- read_csv("final_all_counts.csv") %>%
  select(1, 2:4, 12:14, 19:21)

GeneExampleData <- as.matrix(counts[-1])
row.names(GeneExampleData) = counts$gene

CondVector <- rep(paste("t",1:3,sep=""),each=3)
print(CondVector)



Conditions <- factor(CondVector, levels=c("t1","t2","t3"))
str(Conditions)
levels(Conditions)


Sizes <- MedianNorm(GeneExampleData)

GeneNormData <- GetNormalizedMat(GeneExampleData, Sizes)


PlotExp(GeneNormData, Conditions, Name="ENSMUSG00000025917")

EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=10)

GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)

head(GeneDECalls)
str(GeneDECalls)

GeneConfCalls2 <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.6, OnlyDynamic=TRUE)
str(GeneConfCalls2$EachPathNames)


dd <- unlist(GeneConfCalls2$EachPathNames["Down-Down"])
uu <- unlist(GeneConfCalls2$EachPathNames["Up-Up"])
du <- unlist(GeneConfCalls2$EachPathNames["Down-Up"])
ud <- unlist(GeneConfCalls2$EachPathNames["Up-Down"])



#plotgenes function

norm_counts <- as.data.frame(GeneNormData)
norm_counts$gene <- counts$gene

colnames(norm_counts)[-10] = CondVector

library(janitor)



d1 <- norm_counts %>%
  clean_names() %>%
  pivot_longer(cols = 1:3,
               names_to = "t1",
               values_to = "t1_count") %>%
  select(9,7)

d2 <- norm_counts %>%
  clean_names() %>%
  pivot_longer(cols = 4:6,
               names_to = "t2",
               values_to = "t2_count") %>%
  select(9)


d3 <- norm_counts %>%
  clean_names() %>%
  pivot_longer(cols = 7:9,
               names_to = "t3",
               values_to = "t3_count")  %>%
  select(9)


norm_df <- cbind(d1, d2, d3)

sig_genes <- norm_df %>%
  dplyr::group_by(gene) %>%
  summarise(P30 = median(t1_count),
            P90 = median(t2_count),
            P210 = median(t3_count)) %>%
  mutate(signi = case_when(
    gene %in% dd ~ "Down-Down",
    gene %in% uu ~ "upup",
    gene %in% ud ~ "updown",
    gene %in% du ~ "Down-Up",
    TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(signi)) %>%
  pivot_longer(cols = 2:4,
               names_to = "timepoint",
               values_to = "n_expression")

sig_genes$timepoint <- factor(sig_genes$timepoint, levels = c("P30", "P90", "P210"))


counts <- sig_genes %>%
  dplyr::count(signi) %>%
  mutate(label = as.character(n/3)
  )
# 
# sig_genes <- sig_genes %>%
#   mutate(mygene = case_when(
#     gene == "ENSMUSG00000026678" ~ "RGS5",
#     gene == "ENSMUSG00000030209" ~ "GRIN2B",
#     TRUE ~ NA_character_))

sig_genes <- sig_genes %>%
  mutate(mygene = case_when(
    gene %in% c("ENSMUSG00000051228",
                "ENSMUSG00000000617",
                "ENSMUSG00000030523",
                "ENSMUSG00000070337",
                "ENSMUSG00000021719",
                "ENSMUSG00000024186",
                "ENSMUSG00000025739",
                "ENSMUSG00000004630",
                "ENSMUSG00000031748",
                "ENSMUSG00000023439",
                'ENSMUSG00000056043',
                'ENSMUSG00000032192',
                "ENSMUSG00000015968") ~ "Post-synaptic",
    gene %in% c("ENSMUSG00000093865",
                "ENSMUSG00000048988",
                "ENSMUSG00000043460",
                "ENSMUSG00000039952",
                "ENSMUSG00000042961",
                "ENSMUSG00000031142",
                "ENSMUSG00000045103") ~ "Elfn1",
    gene %in% c("ENSMUSG00000024842,
                ENSMUSG00000041460") ~ "CSNB",
    TRUE ~ NA_character_))

mygenes <- sig_genes %>%
  filter(!is.na(mygene))



(fig1 <- sig_genes %>%
    filter(is.na(mygene)) %>%
    ggplot(aes(x = timepoint,
               y = log(n_expression), group = gene)) +
    theme_pubr() +
    geom_line(alpha = 0.2, color = "grey") +
    geom_point(alpha = 0.2, color = "grey") +
    facet_wrap(vars(signi), ncol =1, drop = TRUE, scales = "free_y") +
    geom_line(data = mygenes, aes(x = timepoint,
                                  y = log(n_expression), group = mygene, color = mygene)) +
    geom_point(data = mygenes, aes(x = timepoint,
                                   y = log(n_expression), group = mygene, color = mygene)) +
    geom_text(data    = counts,
              mapping = aes(x = Inf, y = Inf, label = label),
              hjust   = 1,
              vjust   = 1,
              inherit.aes=FALSE
    ) +  
    geom_text(data    = counts,
              mapping = aes(x = Inf, y = Inf, label = label),
              hjust   = 1,
              vjust   = 1,
              inherit.aes=FALSE
    ) +
    labs(title = "Control",
         y = "Log of Normalized Expression",
         x = "Time Point",
         color = "")) 




(controlfig <- sig_genes %>%
    ggplot(aes(x = timepoint,
               y = log(n_expression), group = gene)) +
    theme_pubr() +
    geom_line(alpha = 0.2, color = "grey") +
    geom_point(alpha = 0.1, color = "grey") +
    facet_wrap(vars(signi), ncol = 1, scales = "free_y") +
    theme(legend.position = "none") +
    geom_text(data    = counts,
              mapping = aes(x = Inf, y = Inf, label = label),
              hjust   = 1,
              vjust   = 1,
              inherit.aes=FALSE
    ) +  
    labs(title = "Control",
         y = "Log of Normalized Expression",
         x = "Time Point"
    ))

library(patchwork)

(fig1 + controlfig )/ guide_area() + plot_annotation(tag_levels = 'A') +  plot_layout(guides = 'collect', nrow = 2)

down <- sig_genes %>%
  filter(signi %in% c("downdown"),
         timepoint %in% c("P30", "P210")) %>%
  pivot_wider(id_cols = gene, names_from  = timepoint, values_from = n_expression) %>%
  mutate(log2_fold_change = log2(P30/P210)) %>%
  arrange(desc(log2_fold_change), desc(P30)) %>%
  slice(1:10) %>%
  mutate(change = "downregulation")

up <- sig_genes %>%
  filter(signi %in% c("upup"),
         timepoint %in% c("P30", "P210")) %>%
  pivot_wider(id_cols = gene, names_from  = timepoint, values_from = n_expression) %>%
  mutate(log2_fold_change = log2(P30/P210)) %>%
  arrange(log2_fold_change, desc(P210)) %>%
  slice(1:10) %>%
  mutate(change = "upregulation")

controls<- rbind(up, down)
controls$condition <- "control"

write.csv(rbind(cng, controls), "ebseq_tcgenes.csv")

