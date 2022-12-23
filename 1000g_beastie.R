library(tidyverse)


#####
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
#####

#load in scarletts pre annotated and filtered data
data_filtered_anno <- read_tsv("/hpc/group/allenlab/scarlett/output/RNAseq/from_little_damon/445_samples_output_cleaned_geneset_1total_2ind_5anc.tsv") %>% 
  filter(totalCount >= 30)

#####
# data_comp_length <- data_filtered_anno %>% mutate(reads_bin = ifelse(totalCount >= 30, T, F))
# 
# ggplot(data_comp_length %>% filter(RR_region != "normal"), 
#        aes(x = RR_region , y = abslog2_posterior_mean, 
#            color = reads_bin)) + geom_boxplot()

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

#####
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
#####
#geneset annotations
geneset_mat <- matrix(data_filtered_anno %>% 
                        select(contains("geneset")) %>% colnames, 
                      ncol = 2, byrow = T)
library(nonpar)

#####
#in gene set vs out of gene set
sum_stats_1 <- tibble(vars = as.character(), facs = as.character(), 
                      mean_in = as.numeric(), 
                      mean_out = as.numeric(), 
                      median_in = as.numeric(),
                      median_out = as.numeric(), 
                      var_in = as.numeric(), 
                      var_out = as.numeric(),
                      sd_in = as.numeric(),
                      sd_out = as.numeric(),
                      .rows = 0)

for(i in 1:nrow(geneset_mat)){
  var <- geneset_mat[i,1]
  fac <- data_filtered_anno %>% filter(get(var) != "other genes")  %>% .[, var] %>% unique
  show(ggplot(data_filtered_anno %>% mutate(across(var, factor, levels = c(fac, "other genes"))),
              aes(x = get(var), y = abslog2_posterior_mean)) +
    geom_boxplot() + xlab(var) + theme_bw() + 
      ggpubr::stat_compare_means(method = "wilcox.test", method.args = list(exact = T)) +
    labs(title = var))
  
  # STAT <- data_filtered_anno  %>% group_by(get(var)) %>% 
  #   summarize(Avg = mean(abslog2_posterior_mean), Median = median(abslog2_posterior_mean)) %>% 
  #   pivot_longer(Avg:Median, names_to = "Stat", values_to = "Value")
  # 
  show(ggplot(data_filtered_anno  %>% mutate(across(var, factor, levels = c(fac, "other genes"))),
              aes(x = abslog2_posterior_mean)) +
         geom_histogram(bins = 100) + theme_bw() +
         #geom_vline(data = STAT, mapping = aes(xintercept = Value, color = Stat)) +
         facet_wrap(~get(var), ncol=2) + labs(title = var))

  mean_1 <- mean(data_filtered_anno %>% filter(get(var)== fac%>%as.character()) %>% .$abslog2_posterior_mean)
  mean_2 <- mean(data_filtered_anno %>% filter(get(var)== "other genes") %>% .$abslog2_posterior_mean)

  median_1 <- median(data_filtered_anno %>% filter(get(var)== fac%>%as.character()) %>% .$abslog2_posterior_mean)
  median_2 <- median(data_filtered_anno %>% filter(get(var)== "other genes") %>% .$abslog2_posterior_mean)
  
  var_1 <- var(data_filtered_anno %>% filter(get(var)== fac%>%as.character()) %>% .$abslog2_posterior_mean)
  var_2 <- var(data_filtered_anno %>% filter(get(var)== "other genes") %>% .$abslog2_posterior_mean)
  
  sd_1 <- sd(data_filtered_anno %>% filter(get(var)== fac%>%as.character()) %>% .$abslog2_posterior_mean)
  sd_2 <- sd(data_filtered_anno %>% filter(get(var)== "other genes") %>% .$abslog2_posterior_mean)

  lepage.test(data_filtered_anno %>% filter(get(var)== fac%>%as.character()) %>% .$abslog2_posterior_mean %>% .[1:50],
              data_filtered_anno %>% filter(get(var)== "other genes") %>% .$abslog2_posterior_mean %>% .[1:50],
              method = "Exact")

  cucconi.test(data_filtered_anno %>% filter(get(var)== fac%>%as.character()) %>% .$abslog2_posterior_mean,
               data_filtered_anno %>% filter(get(var)== "other genes") %>% .$abslog2_posterior_mean,
               method = "permutation")
  
  sum_stats_1 <- add_row(.data = sum_stats_1, vars = var, facs = fac %>% as.character,
                        mean_in = mean_1, mean_out = mean_2,
                        median_in = median_1, median_out = median_2,
                        var_in = var_1, var_out = var_2,
                        sd_in = sd_1, sd_out = sd_2)
}

write_tsv(sum_stats_1, "~/ase_models/sum_stats_1.tsv")

#####
#within gene set, comparing across RR regions

shift_test <- tibble(titles = as.character(), vars = as.character(), 
                     lepage_normal_HRR = as.numeric(), 
                     lepage_normal_CS = as.numeric(), 
                     lepage_CS_HRR = as.numeric(),
                     cucconi_normal_HRR = as.numeric(), 
                     cucconi_normal_CS = as.numeric(), 
                     cucconi_CS_HRR = as.numeric(),
                     .rows = 0)

sum_stats_2 <- tibble(titles = as.character(), vars = as.character(), 
                      mean_normal = as.numeric(), 
                      mean_HRR = as.numeric(), 
                      mean_CS = as.numeric(), 
                      median_normal = as.numeric(), 
                      median_HRR = as.numeric(), 
                      median_CS = as.numeric(),
                      var_normal = as.numeric(), 
                      var_HRR = as.numeric(), 
                      var_CS = as.numeric(),
                      sd_normal = as.numeric(), 
                      sd_HRR = as.numeric(), 
                      sd_CS = as.numeric(),
                      .rows = 0)

for(i in 1:nrow(geneset_mat)){
  title_ <- geneset_mat[i,1]
  var <- geneset_mat[i,2]
  this_gene_set <- data_filtered_anno %>% filter(get(var) == 1)
  # show(ggplot(this_gene_set, #%>% filter(RR_region != "normal"),
  #        aes(x = RR_region, y = abslog2_posterior_mean)) + geom_boxplot() +
  #   theme_bw() + labs(title = title_)) #+ ggpubr::stat_compare_means(method = "wilcox.test", method.args = list(exact = T)))
  
  print(title_)
  print(ggpubr::compare_means(data = this_gene_set, 
                              formula = abslog2_posterior_mean ~ RR_region,
                              method = "wilcox.test",
                              exact = T))
  print("-----")
  # 
  # show(ggplot(this_gene_set, aes(x = abslog2_posterior_mean)) + 
  #   geom_histogram(bins = 100) + theme_bw() + facet_wrap(~RR_region, nrow = 3) + 
  #   labs(title = title_))
  
  # p1 <- try(lepage.test(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean,
  #                       this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean) %>% .$p.value)
  # 
  # p2 <- try(lepage.test(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean,
  #                       this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean) %>% .$p.value)
  # 
  # p3 <- try(lepage.test(this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean,
  #                       this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean) %>% .$p.value)
  # 
  # p4 <- try(cucconi.test(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean,
  #                        this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean) %>% .$p.value)
  # 
  # p5 <- try(cucconi.test(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean,
  #                        this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean) %>% .$p.value)
  # 
  # p6 <- try(cucconi.test(this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean,
  #                        this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean) %>% .$p.value)
  # 
  # shift_test <- add_row(.data = shift_test, titles = title_, vars = var, 
  #                       lepage_normal_HRR = ifelse(class(p1) == "try-error", NA, p1), 
  #                       lepage_normal_CS = ifelse(class(p2) == "try-error", NA, p2), 
  #                       lepage_CS_HRR = ifelse(class(p3) == "try-error", NA, p3), 
  #                       cucconi_normal_HRR = ifelse(class(p4) == "try-error", NA, p4), 
  #                       cucconi_normal_CS = ifelse(class(p5) == "try-error", NA, p5), 
  #                       cucconi_CS_HRR = ifelse(class(p6) == "try-error", NA, p6))
  
  # mean_1 <- mean(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean)
  # mean_2 <- mean(this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean)
  # mean_3 <- mean(this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean)
  # 
  # median_1 <- median(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean)
  # median_2 <- median(this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean)
  # median_3 <- median(this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean)
  # 
  # var_1 <- var(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean)
  # var_2 <- var(this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean)
  # var_3 <- var(this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean)
  # 
  # sd_1 <- sd(this_gene_set %>% filter(RR_region == "normal") %>% .$abslog2_posterior_mean)
  # sd_2 <- sd(this_gene_set %>% filter(RR_region == "HRR") %>% .$abslog2_posterior_mean)
  # sd_3 <- sd(this_gene_set %>% filter(RR_region == "CS") %>% .$abslog2_posterior_mean)
  # 
  # sum_stats_2 <- add_row(.data = sum_stats_2, titles = title_, vars = var,
  #                        mean_normal = mean_1, 
  #                        mean_HRR = mean_2, 
  #                        mean_CS = mean_3, 
  #                        median_normal = median_1, 
  #                        median_HRR = median_2, 
  #                        median_CS = median_3,
  #                        var_normal = var_1, 
  #                        var_HRR = var_2, 
  #                        var_CS = var_3,
  #                        sd_normal = sd_1, 
  #                        sd_HRR = sd_2, 
  #                        sd_CS = sd_3)
  
}

write_tsv(shift_test, "~/ase_models/shift_test.tsv")
write_tsv(sum_stats_2, "~/ase_models/sum_stats_2.tsv")

#####
#what are the genes in the imprinted gene set, specifically those with high ASE in CS and HRR?
#first imprinted gene set: 
  #40 genes in total
data_filtered_anno %>% filter(geneset3 == 1) %>% select(Gene.name) %>% distinct()
  #25 genes with ASE
data_filtered_anno %>% filter(geneset3 == 1, ASE_genes == "ASE") %>% select(Gene.name) %>% unique
  #or 28 using the abslog2 metric
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 1) %>% select(Gene.name) %>% unique
  #9 in CS
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 1, RR_region == "CS") %>% select(Gene.name) %>% unique
  #8 in HRR
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 1, RR_region == "HRR") %>% select(Gene.name) %>% unique
  #11 in Normal
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 1, RR_region == "normal") %>% select(Gene.name) %>% unique
#using higher threshold:
  #5 in CS
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 2, RR_region == "CS") %>% select(Gene.name) %>% unique
  #2 in HRR
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 2, RR_region == "HRR") %>% select(Gene.name) %>% unique
  #3 in Normal
data_filtered_anno %>% filter(geneset3 == 1, abslog2_posterior_mean > 2, RR_region == "normal") %>% select(Gene.name) %>% unique

this_gene_set <- data_filtered_anno %>% filter(geneset3 == 1)
ggplot(this_gene_set, aes(x = abslog2_posterior_mean)) +
  geom_histogram(bins = 100) + theme_bw() + facet_wrap(~RR_region, nrow = 3)
ggplot(this_gene_set %>% filter(abslog2_posterior_mean > 2), 
       aes(x = abslog2_posterior_mean, fill = Gene.name)) +
  geom_histogram(bins = 100) + theme_bw() + facet_wrap(~RR_region, nrow = 3)


#second imprinted set
  #2 genes in total
data_filtered_anno %>% filter(geneset8 == 1) %>% select(Gene.name) %>% distinct()
  #2 genes with ASE
data_filtered_anno %>% filter(geneset8 == 1, ASE_genes == "ASE") %>% select(Gene.name) %>% unique
  #or 2 using the abslog2 metric
data_filtered_anno %>% filter(geneset8 == 1, abslog2_posterior_mean > 1) %>% select(Gene.name) %>% unique
this_gene_set <- data_filtered_anno %>% filter(geneset8 == 1)
ggplot(this_gene_set, aes(x = abslog2_posterior_mean)) +
  geom_histogram(bins = 100) + theme_bw() + facet_wrap(~RR_region, nrow = 3)
ggplot(this_gene_set, 
       aes(x = abslog2_posterior_mean, fill = Gene.name)) +
  geom_histogram(bins = 100) + theme_bw() + facet_wrap(~RR_region, nrow = 3)

#####
### pli and missense scores
library(readxl)
any_ase <- data_filtered_anno %>% filter(ASE_genes == "ASE") %>% .$geneID %>% unique
gene_level_ase <- data_filtered_anno %>% 
  mutate(has_ASE = ifelse(geneID %in% any_ase, T, F)) %>% 
  select(geneID, Gene.stable.ID, Gene.name, has_ASE) %>% distinct()

pli_scores <- read_excel("/hpc/group/allenlab/gene_sets/genes_EDS_RVIS_PLI.xlsx")
hist(pli_scores$pLI %>% as.numeric(), xlab = "pLI", main = "pLI score frequency")

pli_ase <- left_join(gene_level_ase, pli_scores, by = c("Gene.stable.ID" = "GeneSymbol"))
pli_ase$pLI <- pli_ase$pLI %>% as.numeric()
ggplot(pli_ase, aes(x = has_ASE, y = pLI)) + geom_boxplot() + theme_bw() +
  ggpubr::stat_compare_means()

gnomad <- read_tsv("ase_models/supplementary_dataset_11_full_constraint_metrics.tsv.gz") %>% filter(canonical == TRUE)

gnomad_ase <- left_join(gene_level_ase, 
                        gnomad %>% select(gene, transcript, pLI, mis_z), by = c("Gene.name" = "gene"))

hist(gnomad_ase$pLI %>% as.numeric(), xlab = "pLI", main = "pLI score frequency")

ggplot(gnomad_ase, aes(x = has_ASE, y = pLI)) + geom_boxplot() + theme_bw() +
  ggpubr::stat_compare_means()
 
#mis sense z
hist(gnomad_ase$mis_z %>% as.numeric(), xlab = "mis Z", main = "mis Z score frequency")

ggplot(gnomad_ase, aes(x = has_ASE, y = mis_z)) + geom_boxplot() + theme_bw() +
  ggpubr::stat_compare_means()


#haploinsufficient gene sets:
data_filtered_anno %>% filter(geneset1 == 1 | geneset2 == 1)
#####
#look at top genes for all gene sets (ie are just a couple genes contributing to a shift in ASE)

for(i in 1:nrow(geneset_mat)){
  title_ <- geneset_mat[i,1]
  var <- geneset_mat[i,2]
  this_gene_set <- data_filtered_anno %>% filter(get(var) == 1)
  
  print(title_)
  print(var)
  print(paste0("Number of genes in geneset in total: ", 
               this_gene_set %>% select(Gene.name) %>% distinct() %>% nrow()))
  print(paste0("Number of genes in geneset with ASE: ", 
               this_gene_set %>% filter(ASE_genes == "ASE") %>% select(Gene.name) %>% distinct() %>% nrow(),
               " using binary calls"))
  print(paste0("Number of genes in geneset with ASE: ", 
               this_gene_set %>% filter(abslog2_posterior_mean > 3) %>% select(Gene.name) %>% distinct() %>% nrow(),
               " using the abslog2 metric"))
  print(paste0("Number of ASE genes in geneset in CS: ", 
               this_gene_set %>% filter(abslog2_posterior_mean > 3, RR_region == "CS") %>% 
                 select(Gene.name) %>% distinct() %>% nrow()))
  print(paste0("Number of ASE genes in geneset in HRR: ", 
               this_gene_set %>% filter(abslog2_posterior_mean > 3, RR_region == "HRR") %>% 
                 select(Gene.name) %>% distinct() %>% nrow()))
  print(paste0("Number of ASE genes in geneset in Normal: ", 
               this_gene_set %>% filter(abslog2_posterior_mean > 3, RR_region == "normal") %>% 
                 select(Gene.name) %>% distinct() %>% nrow()))
  
  show(ggplot(this_gene_set, aes(x = abslog2_posterior_mean)) +
         geom_histogram(bins = 100) + theme_bw() + 
         facet_wrap(~RR_region, nrow = 3) + labs(title = title_))
  
  show(ggplot(this_gene_set %>% filter(abslog2_posterior_mean > 3), 
              aes(x = abslog2_posterior_mean, fill = Gene.name)) +
         geom_histogram(bins = 100) + theme_bw() + 
         facet_wrap(~RR_region, nrow = 3) + labs(title = title_))
  
}
#####
## gene level summaries
gene_level <- data_filtered_anno %>% group_by(geneID) %>%
  summarise(gene_mean_ase = mean(abslog2_posterior_mean), 
            gene_median_ase = median(abslog2_posterior_mean),
            gene_var_ase = var(abslog2_posterior_mean))

gene_level_sets <- left_join(gene_level, data_filtered_anno %>% select(1, RR_region, contains("geneset")) %>% distinct())

#####
sum_stats_1 <- tibble(vars = as.character(), facs = as.character(), 
                      mean_in = as.numeric(), 
                      mean_out = as.numeric(), 
                      median_in = as.numeric(),
                      median_out = as.numeric(), 
                      var_in = as.numeric(), 
                      var_out = as.numeric(),
                      sd_in = as.numeric(),
                      sd_out = as.numeric(),
                      c_test = as.numeric(),
                      .rows = 0)

for(i in 1:nrow(geneset_mat)){
  var <- geneset_mat[i,1]
  fac <- gene_level_sets %>% filter(get(var) != "other genes")  %>% .[, var] %>% unique
  show(ggplot(gene_level_sets %>% mutate(across(var, factor, levels = c(fac, "other genes"))),
              aes(x = get(var), y = gene_var_ase)) +
         geom_boxplot() + xlab(var) + theme_bw() + 
         ggpubr::stat_compare_means() +
         #ggpubr::stat_compare_means(method = "wilcox.test", method.args = list(exact = T)) +
         labs(title = var))
  
  # show(ggplot(gene_level_sets %>% mutate(across(var, factor, levels = c(fac, "other genes"))),
  #                         aes(x = gene_median_ase)) +
  #                    geom_histogram(bins = 100) + theme_bw() +
  #                    facet_wrap(~get(var), nrow=2) + labs(title = var))
  
  mean_1 <- mean(gene_level_sets %>% filter(get(var)== fac%>%as.character()) %>% .$gene_var_ase)
  mean_2 <- mean(gene_level_sets %>% filter(get(var)== "other genes") %>% .$gene_var_ase)
  
  median_1 <- median(gene_level_sets %>% filter(get(var)== fac%>%as.character()) %>% .$gene_var_ase)
  median_2 <- median(gene_level_sets %>% filter(get(var)== "other genes") %>% .$gene_var_ase)
  
  var_1 <- var(gene_level_sets %>% filter(get(var)== fac%>%as.character()) %>% .$gene_var_ase)
  var_2 <- var(gene_level_sets %>% filter(get(var)== "other genes") %>% .$gene_var_ase)
  
  sd_1 <- sd(gene_level_sets %>% filter(get(var)== fac%>%as.character()) %>% .$gene_var_ase)
  sd_2 <- sd(gene_level_sets %>% filter(get(var)== "other genes") %>% .$gene_var_ase)
  
  # lepage.test(gene_level_sets %>% filter(get(var)== fac%>%as.character()) %>% .$gene_median_ase,
  #             gene_level_sets %>% filter(get(var)== "other genes") %>% .$gene_median_ase,
  #             method = "Exact")
  
  c_test <- cucconi.test(gene_level_sets %>% filter(get(var)== fac%>%as.character()) %>% .$gene_var_ase,
                         gene_level_sets %>% filter(get(var)== "other genes") %>% .$gene_var_ase,
                         method = "permutation")
  
  sum_stats_1 <- add_row(.data = sum_stats_1, vars = var, facs = fac %>% as.character,
                         mean_in = mean_1, mean_out = mean_2,
                         median_in = median_1, median_out = median_2,
                         var_in = var_1, var_out = var_2,
                         sd_in = sd_1, sd_out = sd_2, c_test = c_test$p.value)
}

write_tsv(sum_stats_1, "ase_models/sum_stats_1_gene_variance.tsv")

#####
shift_test <- tibble(titles = as.character(), vars = as.character(), 
                     lepage_normal_HRR = as.numeric(), 
                     lepage_normal_CS = as.numeric(), 
                     lepage_CS_HRR = as.numeric(),
                     cucconi_normal_HRR = as.numeric(), 
                     cucconi_normal_CS = as.numeric(), 
                     cucconi_CS_HRR = as.numeric(),
                     .rows = 0)

sum_stats_2 <- tibble(titles = as.character(), vars = as.character(), 
                      mean_normal = as.numeric(), 
                      mean_HRR = as.numeric(), 
                      mean_CS = as.numeric(), 
                      median_normal = as.numeric(), 
                      median_HRR = as.numeric(), 
                      median_CS = as.numeric(),
                      var_normal = as.numeric(), 
                      var_HRR = as.numeric(), 
                      var_CS = as.numeric(),
                      sd_normal = as.numeric(), 
                      sd_HRR = as.numeric(), 
                      sd_CS = as.numeric(),
                      .rows = 0)

for(i in 1:nrow(geneset_mat)){
  title_ <- geneset_mat[i,1]
  var <- geneset_mat[i,2]
  this_gene_set <- gene_level_sets %>% filter(get(var) == 1)
  # show(ggplot(this_gene_set, #%>% filter(RR_region != "normal"),
  #        aes(x = RR_region, y = gene_median_ase)) + geom_boxplot() +
  #   theme_bw() + labs(title = title_) + ggpubr::stat_compare_means())

  print(title_)
  print(ggpubr::compare_means(data = this_gene_set,
                              formula = gene_median_ase ~ RR_region,
                              method = "wilcox.test",
                              exact = T))
  print("-----")
  
  
  # show(ggplot(this_gene_set, aes(x = gene_median_ase)) +
  #   geom_histogram(bins = 100) + theme_bw() + facet_wrap(~RR_region, nrow = 3) +
  #   labs(title = title_))
  # 
  # p1 <- try(lepage.test(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase,
  #                       this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase) %>% .$p.value)
  # 
  # p2 <- try(lepage.test(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase,
  #                       this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase) %>% .$p.value)
  # 
  # p3 <- try(lepage.test(this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase,
  #                       this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase) %>% .$p.value)
  # 
  # p4 <- try(cucconi.test(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase,
  #                        this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase) %>% .$p.value)
  # 
  # p5 <- try(cucconi.test(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase,
  #                        this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase) %>% .$p.value)
  # 
  # p6 <- try(cucconi.test(this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase,
  #                        this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase) %>% .$p.value)
  # 
  # shift_test <- add_row(.data = shift_test, titles = title_, vars = var,
  #                       lepage_normal_HRR = ifelse(class(p1) == "try-error", NA, p1),
  #                       lepage_normal_CS = ifelse(class(p2) == "try-error", NA, p2),
  #                       lepage_CS_HRR = ifelse(class(p3) == "try-error", NA, p3),
  #                       cucconi_normal_HRR = ifelse(class(p4) == "try-error", NA, p4),
  #                       cucconi_normal_CS = ifelse(class(p5) == "try-error", NA, p5),
  #                       cucconi_CS_HRR = ifelse(class(p6) == "try-error", NA, p6))
  # 
  # mean_1 <- mean(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase)
  # mean_2 <- mean(this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase)
  # mean_3 <- mean(this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase)
  # 
  # median_1 <- median(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase)
  # median_2 <- median(this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase)
  # median_3 <- median(this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase)
  # 
  # var_1 <- var(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase)
  # var_2 <- var(this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase)
  # var_3 <- var(this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase)
  # 
  # sd_1 <- sd(this_gene_set %>% filter(RR_region == "normal") %>% .$gene_median_ase)
  # sd_2 <- sd(this_gene_set %>% filter(RR_region == "HRR") %>% .$gene_median_ase)
  # sd_3 <- sd(this_gene_set %>% filter(RR_region == "CS") %>% .$gene_median_ase)
  # 
  # sum_stats_2 <- add_row(.data = sum_stats_2, titles = title_, vars = var,
  #                        mean_normal = mean_1,
  #                        mean_HRR = mean_2,
  #                        mean_CS = mean_3,
  #                        median_normal = median_1,
  #                        median_HRR = median_2,
  #                        median_CS = median_3,
  #                        var_normal = var_1,
  #                        var_HRR = var_2,
  #                        var_CS = var_3,
  #                        sd_normal = sd_1,
  #                        sd_HRR = sd_2,
  #                        sd_CS = sd_3)
  
}

# write_tsv(shift_test, "~/ase_models/shift_test_gene_median.tsv")
# write_tsv(sum_stats_2, "~/ase_models/sum_stats_2_gene_median.tsv")


#####
##mixed linear model

lm_data <- data_filtered_anno %>% select(sample, 
                                         geneID,
                                         abslog2_posterior_mean, 
                                         RR_region, contains("geneset")) %>% 
  mutate(sample_num = str_sub(sample, start = 3) %>% as.numeric(),
         RR_region = factor(x = RR_region, levels = c("normal", "CS", "HRR")))

library(lme4)
geneset_cols <- paste0("geneset", c(seq(1,8), 10))
ase_lm_null = lmer(log(abslog2_posterior_mean) ~ (1 | geneID), data = lm_data)

#geneset models
ase_lm_list <- list()
for(i in geneset_cols){
  ase_lm = lmer(log(abslog2_posterior_mean) ~ get(i) + (1 | geneID), data = lm_data)
  ase_lm_list[[i]] <- ase_lm
}

p_values <- tibble(geneset = as.character(), 
                   intercept = as.numeric(), 
                   beta = as.numeric(), 
                   pvalue = as.numeric(), .rows = 0)

for(i in geneset_cols){
  intercept <- summary(ase_lm_list[[i]])$coefficients[1]
  beta <-  summary(ase_lm_list[[i]])$coefficients[2]
  anova <- anova(ase_lm_null, ase_lm_list[[i]])
  p_values <- add_row(.data = p_values, 
                      geneset = i, 
                      intercept = intercept,
                      beta = beta,
                      pvalue = anova$`Pr(>Chisq)`[2])
}

#rr region model
ase_lm_rr = lmer(log(abslog2_posterior_mean) ~ RR_region + (1 | geneID), data = lm_data)
anova(ase_lm_null, ase_lm_rr)

#interaction models
ase_lm_interaction_null_list <- list()
ase_lm_interaction_list <- list()
for(i in geneset_cols){
  ase_lm_null_int = lmer(log(abslog2_posterior_mean) ~ get(i) + RR_region + (1 | geneID), data = lm_data)
  ase_lm = lmer(log(abslog2_posterior_mean) ~ get(i) + RR_region + get(i) * RR_region + (1 | geneID), data = lm_data)
  
  ase_lm_interaction_null_list[[i]] <- ase_lm_null_int
  ase_lm_interaction_list[[i]] <- ase_lm
}

p_values <- tibble(geneset = as.character(), 
                   null_intercept = as.numeric(), 
                   null_geneset_beta = as.numeric(), 
                   null_RR_CS_beta = as.numeric(), 
                   null_RR_HRR_beta = as.numeric(), 
                   model_intercept = as.numeric(), 
                   model_geneset_beta = as.numeric(), 
                   model_RR_CS_beta = as.numeric(), 
                   model_RR_HRR_beta = as.numeric(),
                   geneset_CS_interaction = as.numeric(),
                   geneset_HRR_interaction = as.numeric(),
                   pvalue = as.numeric(), .rows = 0)

for(i in geneset_cols){
  null_intercept <- summary(ase_lm_interaction_null_list[[i]])$coefficients[1]
  null_geneset_beta <- summary(ase_lm_interaction_null_list[[i]])$coefficients[2]
  null_RR_CS_beta <- summary(ase_lm_interaction_null_list[[i]])$coefficients[3]
  null_RR_HRR_beta <- summary(ase_lm_interaction_null_list[[i]])$coefficients[4]
  
  model_intercept <- summary(ase_lm_interaction_list[[i]])$coefficients[1]
  model_geneset_beta <- summary(ase_lm_interaction_list[[i]])$coefficients[2]
  model_RR_CS_beta <- summary(ase_lm_interaction_list[[i]])$coefficients[3]
  model_RR_HRR_beta <- summary(ase_lm_interaction_list[[i]])$coefficients[4]
  
  geneset_CS_interaction <- summary(ase_lm_interaction_list[[i]])$coefficients[5]
  geneset_HRR_interaction <- summary(ase_lm_interaction_list[[i]])$coefficients[6]
  
  anova <- anova(ase_lm_interaction_null_list[[i]], ase_lm_interaction_list[[i]])
  
  p_values <- add_row(.data = p_values, 
                      geneset = i, 
                      null_intercept = null_intercept, 
                      null_geneset_beta = null_geneset_beta, 
                      null_RR_CS_beta = null_RR_CS_beta, 
                      null_RR_HRR_beta = null_RR_HRR_beta, 
                      model_intercept = model_intercept, 
                      model_geneset_beta = model_geneset_beta, 
                      model_RR_CS_beta = model_RR_CS_beta, 
                      model_RR_HRR_beta = model_RR_HRR_beta,
                      geneset_CS_interaction = geneset_CS_interaction,
                      geneset_HRR_interaction = geneset_HRR_interaction,
                      pvalue = anova$`Pr(>Chisq)`[2])
}

