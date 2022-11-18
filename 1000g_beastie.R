library(tidyverse)

#this code loads in scarlett's output and does filtering and annotating. 
#code below uses the data scarlett already filtered and annotated.

# scarlett_data <- read_tsv("/hpc/group/allenlab/scarlett/output/RNAseq/445_samples_output.tsv") %>%
#   mutate(ensembl_gene_id = str_split(string = geneID, pattern = "\\.", simplify = T) %>% .[,1])
# #mapping_data <- read_tsv("/hpc/group/allenlab/scarlett/recombination_map_Supp/CS_HRR_exome_final.bed")
# recomb_data <- read_tsv("/hpc/group/allenlab/scarlett/recombination_map_Supp/CS_exome_rates.txt")
# snv_data <- read_tsv("/hpc/group/allenlab/scarlett/recombination_map_Supp/FCQ_snvs.annot.gz")
# 
# ##add gene name annotation
# library(biomaRt)
# mart <- useMart("ensembl")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
# filters <- c("ensembl_gene_id")
# values <- list(ensembl_gene_id = scarlett_data$ensembl_gene_id %>% unique)
# attributes <- c("ensembl_gene_id", "external_gene_name",
#                 "chromosome_name", "start_position", "end_position")
# gene_names <- getBM(attributes = attributes, filters = filters,
#                     values = values, mart = mart) %>% as_tibble()
# detach("package:biomaRt", unload = T)
# 
# scarlett_data <- full_join(scarlett_data, gene_names)
# 
# ## filter genes
# #Only keep genes that exist in at least 2 individuals from each ancestry group.
# #Filter reads for gene coverage at 30. (post-step)
# count_genes <- scarlett_data %>% count(geneID, ancestry) %>%
#   pivot_wider(names_from = ancestry, values_from = n)
# keep_genes <- count_genes %>%
#   filter(if_all(everything(), ~ . >= 2)) %>% .$geneID
# 
# data_filtered <- scarlett_data %>% filter(geneID %in% keep_genes,
#                                           totalCount >= 30) %>%
#   select(geneID, external_gene_name, chromosome_name, start_position, end_position,
#          abslog2_posterior_mean, posterior_mass_support_ALT, ASE,
#          sex, ancestry, sample)
# 
# ##overlap gene_names (and location info) with recomb data
# gene_recomb <- tibble(ensembl_gene_id = c(), external_gene_name = c(),
#                       chromosome_name = c(), start_position = c(), end_position = c(),
#                       chr = c(), start = c(), end = c(), name = c(),
#                       rrates_ceu = c(), rrates_yri = c(), rrates_mtl = c(),
#                       decode = c(), .rows = 0)
# #this finds genes that are completely within a recombination region, or partial overlaps (i think)
# for(i in 1:nrow(gene_names)){
#   chr_name <- paste0("chr", gene_names[i,3])
#   start_p <- gene_names[i,4] %>% as.numeric()
#   end_p <- gene_names[i,5] %>% as.numeric()
#   find_recomb <- recomb_data %>% filter(chr == chr_name, start <= end_p, end >= start_p)
#   new_row <- bind_cols(gene_names[i,], find_recomb)
#   gene_recomb <- bind_rows(gene_recomb, new_row)
# }
# 
# data_recomb <- left_join(data_filtered, gene_recomb %>% select(-ensembl_gene_id))


#load in scarletts pre annotated and filtered data
data_filtered_anno <- read_tsv("/hpc/group/allenlab/scarlett/output/RNAseq/from_little_damon/445_samples_output_cleaned_geneset_1total_2ind_5anc.tsv") %>% 
  filter(totalCount >= 30)

data_comp_length <- data_filtered_anno %>% mutate(reads_bin = ifelse(totalCount >= 30, T, F))

ggplot(data_comp_length, aes(x = RR_region != "normal", 
                             y = abslog2_posterior_mean, 
                             color = reads_bin)) + geom_boxplot()

### Do associations

# ggplot(data_recomb %>% filter(ancestry == "CEU"), aes(x = ASE, y = rrates_ceu)) +
#   geom_boxplot() + ggpubr::stat_compare_means() + theme_bw()
# 
# ggplot(data_recomb %>% filter(ancestry == "YRI"), aes(x = ASE, y = rrates_yri)) +
#   geom_boxplot() + ggpubr::stat_compare_means() + theme_bw()
# 
# 
# ggplot(data_recomb %>% filter(ancestry == "CEU"), aes(x = abslog2_posterior_mean,
#                                                       y = rrates_ceu)) +
#   geom_point()
# 
# ggplot(data_recomb %>% filter(ancestry == "YRI"), aes(x = abslog2_posterior_mean,
#                                                       y = rrates_yri)) +
#   geom_point()

ggplot(data_filtered_anno, aes(x = RR_region , y = abslog2_posterior_mean)) +
  geom_boxplot() + theme_bw() + ggpubr::stat_compare_means()

ggplot(data_filtered_anno %>% filter(RR_region != "normal"), 
       aes(x = RR_region , y = abslog2_posterior_mean)) + 
  geom_boxplot() + theme_bw() + ggpubr::stat_compare_means()

ggplot(data_filtered_anno %>% filter(RR_region != "normal"), 
       aes(x = RR_region , y = abslog2_posterior_mean, color = ancestry)) + 
  geom_boxplot() + theme_bw() + ggpubr::stat_compare_means()

gene_level <- data_filtered_anno %>% filter(RR_region != "normal") %>% group_by(geneID, RR_region) %>%
  summarise(gene_mean_ase = mean(abslog2_posterior_mean))

ggplot(gene_level,# %>% filter(gene_mean_ase < 3), 
       aes(x = RR_region , y = gene_mean_ase)) + 
  geom_boxplot() + theme_bw() + ggpubr::stat_compare_means()

gene_level <- data_filtered_anno %>% filter(RR_region != "normal") %>% group_by(geneID, RR_region, ancestry) %>%
  summarise(gene_mean_ase = mean(abslog2_posterior_mean))

ggplot(gene_level, 
       aes(x = RR_region , y = gene_mean_ase, color = ancestry)) + 
  geom_boxplot() + theme_bw() + ggpubr::stat_compare_means()

ggplot(data_filtered_anno %>% filter(RR_region != "normal", geneID != "ENSG00000070756.9"), 
       aes(x = RR_region , y = abslog2_posterior_mean)) + 
  geom_boxplot() + theme_bw() + ggpubr::stat_compare_means()

ggplot(data_filtered_anno %>% filter(RR_region == "CS", geneID != "ENSG00000070756.9", abslog2_posterior_mean > 5), 
       aes(x = RR_region , y = abslog2_posterior_mean, color = geneID)) + 
  geom_boxplot() + theme_bw() #+ ggpubr::stat_compare_means()


data_filtered_anno %>% filter(RR_region == "CS") %>% arrange(desc(abslog2_posterior_mean))
hist(data_filtered_anno %>% filter(RR_region == "CS", abslog2_posterior_mean > 2) %>% .$abslog2_posterior_mean, 
     main = "Freq of high ASE in CS regions", xlab = "abslog2_posterior_mean > 2")
hist(data_filtered_anno %>% filter(RR_region == "HRR", abslog2_posterior_mean > 2) %>% .$abslog2_posterior_mean,
     main = "Freq of high ASE in HRR regions", xlab = "abslog2_posterior_mean > 2")


#In HRR/CS comparison, look at binary ASE values instead of continuous.
ggplot(data_filtered_anno, aes(x = RR_region, y = ASE)) + geom_point()
#how many samples in each category?
data_filtered_anno %>% count(RR_region)
data_filtered_anno %>% count(ASE)
data_filtered_anno %>% count(RR_region, ASE)

##associate ASE with SNV annotations
#there are different values for different SNPs in the genes... not sure which ones to use
# ase_snv <- left_join(data_filtered, snv_data %>% select(gene, functionConsensus),
#                      by = c("external_gene_name" = "gene"))
