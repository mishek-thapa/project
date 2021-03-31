library(EBSeqHMM)
library(tidyverse)
library(ggpubr)

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

