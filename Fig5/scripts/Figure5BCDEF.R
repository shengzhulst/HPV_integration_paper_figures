# Author: Shaomiao Xia
# Description: Analysis of gene-level integration sites and cumulative distribution plots for Figure 5B, 5C, 5D, 5E, and 5F.

# Load required libraries
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(data.table)

# Set working directory
setwd("HPV_integration_paper_figures/Fig5")

# Read input data
WGS_sites <- read.csv("./data/WGS_combined_gene_annotations.csv")
RNA_filtered <- read.csv("./data/all_combined_genes_annot_filtered.csv")
RNA_suspect <- read.csv("./data/all_combined_genes_annot_filtered_sites_suspected.csv")
RNA_filtered$group <- "unspliced"
RNA_suspect$group <- "spliced"
RNA_sites <- rbind(RNA_filtered, RNA_suspect)

# Standardize annotation labels
WGS_sites$annot <- ifelse(WGS_sites$from == "chipenrich", "enhancer_or_nearest_TSS", WGS_sites$annot)
WGS_sites$annot <- ifelse(WGS_sites$from == "uscs", "inside_genes", WGS_sites$annot)

# Read sample list and standardize names
sample_102 <- readLines("./data/samples_for_DR_compare.txt")
sample_102 <- gsub("[.]", "-", sample_102)

# Subset WGS and RNA sites for selected samples
WGS_102 <- WGS_sites[WGS_sites$SampleID %in% sample_102,]
RNA_102 <- RNA_sites[RNA_sites$SampleID %in% sample_102,]

# Filter gene-level sites and remove unwanted annotations
filter_annot <- c("left_two", "left_one", "right_one", "right_two")
WGS_sites_gene <- WGS_102[!WGS_102$annot %in% filter_annot & !is.na(WGS_102$HmGene),] %>% distinct()
RNA_sites_gene <- RNA_102[!RNA_102$annot %in% filter_annot & !is.na(RNA_102$HmGene),] %>% distinct()

# Pivot to wide format and summarize gene annotations
pivot_and_summarize <- function(df) {
  df_wide <- df %>%
    group_by(SampleID, Chr, Start) %>%
    pivot_wider(names_from = annot,
                values_from = HmGene,
                values_fill = NA,
                values_fn = function(x) paste0(x, collapse = ",")) %>%
    group_by(SampleID, Chr, Start) %>%
    summarize(across(.cols = everything(),
                     .fns = ~ paste(unique(na.omit(.x)), collapse = ","),
                     .names = "{.col}"), .groups = "drop")
  colnames(df_wide) <- gsub("hg38_", "", colnames(df_wide))
  df_wide
}

WGS_sites_gene_wide <- pivot_and_summarize(WGS_sites_gene)
RNA_sites_gene_wide <- pivot_and_summarize(RNA_sites_gene)

# Order gene annotations and extract confident gene
annotation_order <- c("genes_firstexons", "genes_exons", "genes_introns",
                     "genes_5UTRs", "genes_promoters", "genes_1to5kb",
                     "genes_3UTRs", "genes_exonintronboundaries",
                     "genes_intronexonboundaries", "genes_cds", "integration",
                     "enhancer_or_nearest_TSS", "inside_genes")

order_annot_gene <- function(df) {
  df$OrderedAnnotGene <- apply(df[, annotation_order], 1, function(x) paste(na.omit(x[x != ""]), collapse = ","))
  df$ConfidentGene <- sapply(strsplit(df$OrderedAnnotGene, ","), function(x) x[1])
  setDT(df)
  df[, Annot := apply(.SD, 1, function(x) paste(names(x)[x != "" & !is.na(x)], collapse = ",")), .SDcols = annotation_order]
  df
}

WGS_sites_gene_wide <- order_annot_gene(WGS_sites_gene_wide)
RNA_sites_gene_wide <- order_annot_gene(RNA_sites_gene_wide)

# Group sites by 20,000 bp range and summarize
group_and_summarize <- function(df) {
  df <- setDT(df)
  df[, group_key := floor(Start / 20000), by = .(SampleID, Chr, ConfidentGene)]
  # Check if 'group' column exists
  if ("group" %in% colnames(df)) {
    df_combined <- df[, .(Start = mean(Start),
                          HPVType = first(HPVType),
                          HPVStart = first(HPVStart),
                          HPVGene = first(HPVGene),
                          HmStrand = first(HmStrand),
                          HPVStrand = first(HPVStrand),
                          Source = first(Source),
                          Annot = first(Annot),
                          OrderedAnnotGene = first(OrderedAnnotGene),
                          from = first(from),
                          group = first(group)),
                      by = .(SampleID, Chr, ConfidentGene, group_key)]
  } else {
    df_combined <- df[, .(Start = mean(Start),
                          HPVType = first(HPVType),
                          HPVStart = first(HPVStart),
                          HPVGene = first(HPVGene),
                          HmStrand = first(HmStrand),
                          HPVStrand = first(HPVStrand),
                          Source = first(Source),
                          Annot = first(Annot),
                          OrderedAnnotGene = first(OrderedAnnotGene),
                          from = first(from)),
                      by = .(SampleID, Chr, ConfidentGene, group_key)]
  }
  df_combined[, group_key := NULL]
  df_combined
}

wgs_combined <- group_and_summarize(WGS_sites_gene_wide)
rna_combined <- group_and_summarize(RNA_sites_gene_wide)

# Distinct sites for merging
WGS_102_distinct <- WGS_102 %>% distinct(SampleID, Chr, Start, .keep_all = TRUE)
RNA_102_distinct <- RNA_102 %>% distinct(SampleID, Chr, Start, .keep_all = TRUE)
WGS_102_distinct_sub <- WGS_102_distinct[, c("SampleID", "Chr", "Start", "HPVType", "HPVStart")]
RNA_102_distinct_sub <- RNA_102_distinct[, c("SampleID", "Chr", "Start", "HPVType", "HPVStart")]

# Merge RNA and WGS sites for distance calculation
RD_merged_df <- merge(RNA_sites_gene_wide, WGS_sites_gene_wide, by = c("SampleID", "Chr"), all.x = TRUE, allow.cartesian = TRUE)
DR_merged_df <- merge(WGS_sites_gene_wide, RNA_sites_gene_wide, by = c("SampleID", "Chr"), all.x = TRUE, allow.cartesian = TRUE)
RD_merged_df$distance <- abs(RD_merged_df$Start.x - RD_merged_df$Start.y)
DR_merged_df$distance <- abs(DR_merged_df$Start.x - DR_merged_df$Start.y)

# Find nearest WGS site for each RNA site
RNA_sites_to_nearest_WGS_sites <- RD_merged_df %>%
  group_by(SampleID, Chr, Start.x) %>%
  slice_min(order_by = distance, n = 1)

# Prepare distance columns for plotting
RNA_sites_to_nearest_WGS_sites$plot_dist <- ifelse(RNA_sites_to_nearest_WGS_sites$distance <= 1e8, RNA_sites_to_nearest_WGS_sites$distance, 1e8)
RNA_sites_to_nearest_WGS_sites$plot_dist <- ifelse(is.na(RNA_sites_to_nearest_WGS_sites$distance), 1e8, RNA_sites_to_nearest_WGS_sites$distance)
RNA_sites_to_nearest_WGS_sites$log_distance <- log10(RNA_sites_to_nearest_WGS_sites$plot_dist + 1)

# Handle overlap group and combine distances
overlap_group <- RNA_sites_to_nearest_WGS_sites[RNA_sites_to_nearest_WGS_sites$group == "unspliced,spliced",]
overlap_group$group <- "spliced"
RNA_sites_to_nearest_WGS_sites$group <- gsub("unspliced,spliced", "unspliced", RNA_sites_to_nearest_WGS_sites$group)
RNA_sites_to_nearest_WGS_sites <- rbind(RNA_sites_to_nearest_WGS_sites, overlap_group)
distance_combined <- RNA_sites_to_nearest_WGS_sites
distance_combined$group <- "combined"
all_groups_distance <- rbind(RNA_sites_to_nearest_WGS_sites, distance_combined)

# Prepare random sites for control
original_site_counts <- data.frame(
  Chr = c(1:22, "X", "Y"),
  Count = c(45, 22, 13, 36, 2, 10, 2, 16, 12, 7, 29, 24, 3, 17, 2, 109, 59, 1, 24, 9, 33, 27, 7, 2)
)
hg38_chr_sizes <- data.frame(
  Chr = c(1:22, "X", "Y"),
  Size = as.numeric(gsub(",", "", c(
    "248,956,422", "242,193,529", "198,295,559", "190,214,555", "181,538,259",
    "170,805,979", "159,345,973", "145,138,636", "138,394,717", "133,797,422",
    "135,086,622", "133,275,309", "114,364,328", "107,043,718", "101,991,189",
    "90,338,345", "83,257,441", "80,373,285", "58,617,616", "64,444,167",
    "46,709,983", "50,818,468", "156,040,895", "57,227,415"
  )))
)

# Generate random sites for all chromosomes
generate_random_sites <- function(chr, count, size) {
  data.frame(Chr = rep(chr, count), Position = sample(1:size, count))
}
random_sites_by_chr <- original_site_counts %>%
  inner_join(hg38_chr_sizes, by = "Chr") %>%
  rowwise() %>%
  do(generate_random_sites(.$Chr, .$Count, .$Size)) %>%
  ungroup()
random_sites_by_chr$Chr <- paste0("chr", random_sites_by_chr$Chr)

# Merge random sites with WGS sites and calculate distances
random_merged_df <- merge(random_sites_by_chr, WGS_102_distinct_sub, by = "Chr", all.x = TRUE)
random_merged_df$distance <- abs(random_merged_df$Start - random_merged_df$Position)
random_sites_to_nearest_WGS_sites <- random_merged_df %>%
  group_by(Chr, Start) %>%
  slice_min(order_by = distance, n = 1)
random_sites_to_nearest_WGS_sites$group <- "random"
random_sites_to_nearest_WGS_sites$plot_dist <- ifelse(random_sites_to_nearest_WGS_sites$distance <= 1e8, random_sites_to_nearest_WGS_sites$distance, 1e8)
random_sites_to_nearest_WGS_sites$plot_dist <- ifelse(is.na(random_sites_to_nearest_WGS_sites$distance), 1e8, random_sites_to_nearest_WGS_sites$distance)
random_sites_to_nearest_WGS_sites$log_distance <- log10(random_sites_to_nearest_WGS_sites$plot_dist + 1)

df_for_cd_plot <- rbind(all_groups_distance[, c("log_distance", "group")], random_sites_to_nearest_WGS_sites[, c("log_distance", "group")])

# Cumulative distribution plot for Figure 5B
svg("./plots/figure5B.svg", width = 5, height = 5)
ggplot(df_for_cd_plot, aes(log_distance, colour = group)) + 
  stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.8) +
  labs(title = " ", x = bquote(log[10](distance+1)), y = "Cumulative percentage of RNA sites (%)", colour = "") +
  scale_x_continuous(limits = c(0, max(df_for_cd_plot$log_distance)), breaks = seq(0, max(df_for_cd_plot$log_distance), by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = function(x) paste0(x * 100)) +
  scale_color_manual(values = c("mediumpurple1", "aquamarine3", "grey70", "indianred1")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = NA),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.1, linetype = "dashed"),
        panel.grid.minor.y = element_line(colour = "grey", linewidth = 0.1, linetype = "dashed"),
        panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA))
dev.off()

# Prepare data for Figure 5C
RD_combined_sub <- all_groups_distance[, c("SampleID", "Chr", "Start.x", "Start.y", "HPVType.x", "HPVType.y", "HPVGene.x", "HPVGene.y", "ConfidentGene.x", "ConfidentGene.y", "OrderedAnnotGene.x", "distance", "log_distance", "group")]
RD_combined_sub$kb_distance <- RD_combined_sub$distance / 1000
RD_combined_sub$RNA_covered_gene <- lapply(strsplit(RD_combined_sub$OrderedAnnotGene.x, ","), function(x) as.vector(unique(x)))
RD_combined_sub$match <- mapply(`%in%`, RD_combined_sub$ConfidentGene.y, RD_combined_sub$RNA_covered_gene)
RD_combined_sub$plot_kb_dist <- ifelse(RD_combined_sub$match == TRUE, RD_combined_sub$kb_distance, 200)
RD_combined_sub$plot_kb_dist <- ifelse(RD_combined_sub$kb_distance <= 200, RD_combined_sub$kb_distance, 200)

# Cumulative distribution plot for Figure 5C
svg("./plots/figure5C.svg", width = 5, height = 5)
ggplot(RD_combined_sub, aes(x = plot_kb_dist, color = group)) +
  stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.8) +
  labs(title = " ", x = bquote(Distance ~ (kb)), y = "Cumulative percentage of matches (%)", colour = "") +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 20)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = function(x) paste0(x * 100)) +
  scale_color_manual(values = c("mediumpurple1", "aquamarine3", "indianred1")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = NA),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.1, linetype = "dashed"),
        panel.grid.minor.y = element_line(colour = "grey", linewidth = 0.1, linetype = "dashed"),
        panel.border = element_rect(linetype = "solid", linewidth = 1, fill = NA))
dev.off()

# PRODUCTIVE SITES ANALYSIS
DR_merged_df <- merge(wgs_combined, rna_combined, by = c("SampleID", "Chr"), all.x = TRUE, allow.cartesian = TRUE)
DR_merged_df$distance <- abs(DR_merged_df$Start.x - DR_merged_df$Start.y)
WGS_sites_to_nearest_RNA_sites <- DR_merged_df %>%
  group_by(SampleID, Chr, Start.x) %>%
  slice_min(order_by = distance, n = 1)
WGS_sites_productive <- WGS_sites_to_nearest_RNA_sites[(WGS_sites_to_nearest_RNA_sites$distance <= 100000) & !is.na(WGS_sites_to_nearest_RNA_sites$distance),]
WGS_sites_silent <- WGS_sites_to_nearest_RNA_sites[WGS_sites_to_nearest_RNA_sites$distance > 100000 | is.na(WGS_sites_to_nearest_RNA_sites$distance),]

prod_silent_count <- data.frame(
  category = c("productive", "silent"),
  count = c(nrow(WGS_sites_productive), nrow(WGS_sites_silent))
)
prod_silent_count$percentage <- prod_silent_count$count / sum(prod_silent_count$count) * 100

# Bar plot for productive/silent sites (Figure 5D)
svg("./plots/figure5D.svg", width = 5, height = 2)
ggplot(prod_silent_count, aes(x = factor(1), y = percentage, fill = category, cumulative = TRUE)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("darkorange1", "skyblue")) +
  geom_text(aes(label = paste0(count, " (", round(percentage, 2), "%)")), position = position_stack(vjust = 0.5)) +
  labs(y = "Percentage", x = "", fill = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.background = element_rect(fill = NA),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", linewidth = 0.8, fill = NA))
dev.off()

# Prepare gene-level summary for productive/silent sites
productive_sites_sub <- WGS_sites_productive[, c("SampleID", "Chr", "Start.x", "HPVType.x", "HPVStart.x", "HPVGene.x", "ConfidentGene.x")]
silent_sites_sub <- WGS_sites_silent[, c("SampleID", "Chr", "Start.x", "HPVType.x", "HPVStart.x", "HPVGene.x", "ConfidentGene.x")]

productive_gene_count <- productive_sites_sub %>%
  group_by(ConfidentGene.x) %>%
  summarize(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))
productive_gene_count$ConfidentGene.x <- factor(productive_gene_count$ConfidentGene.x, levels = productive_gene_count$ConfidentGene.x)

silent_gene_count <- silent_sites_sub %>%
  group_by(ConfidentGene.x) %>%
  summarize(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))
silent_gene_count$ConfidentGene.x <- factor(silent_gene_count$ConfidentGene.x, levels = silent_gene_count$ConfidentGene.x)

# Recurrent loci analysis
recurrent_loci <- read.csv("/Users/shaomiao/Desktop/Umich/SartorLab/HPV_integration/sites/Figure1A_human_recurrent_final_vis_filtered.csv")
loci <- unique(recurrent_loci$HmGene)
productive_sites_sub$group <- "productive"
silent_sites_sub$group <- "silent"
classified_wgs_sites <- rbind(productive_sites_sub, silent_sites_sub)
loci_lookup <- unlist(lapply(strsplit(loci, ":"), function(genes) setNames(rep(paste(genes, collapse = ":"), length(genes)), genes)))
classified_wgs_sites <- classified_wgs_sites %>%
  mutate(loci = ifelse(ConfidentGene.x %in% names(loci_lookup), loci_lookup[ConfidentGene.x], ConfidentGene.x))
classified_wgs_sites$reccurent <- ifelse(classified_wgs_sites$loci %in% loci, TRUE, FALSE)

sites_count <-
  matrix(c(nrow(classified_wgs_sites[classified_wgs_sites$reccurent == TRUE & classified_wgs_sites$group == "productive",]),
           nrow(classified_wgs_sites[classified_wgs_sites$reccurent == FALSE & classified_wgs_sites$group == "productive",]),
           nrow(classified_wgs_sites[classified_wgs_sites$reccurent == TRUE & classified_wgs_sites$group == "silent",]),
           nrow(classified_wgs_sites[classified_wgs_sites$reccurent == FALSE & classified_wgs_sites$group == "silent",])),
         nrow = 2,
         dimnames = list(c("reccurent", "other"),
                         c("productive", "silent")))

fisher.test(sites_count, alternative = "two.sided")

sites_count_m <- reshape2::melt(sites_count, varnames = c("Region", "IntegrationStatus"), id.vars = "IntegrationStatus")

# Bar plot for recurrent/other loci (Figure 5E)
svg("./plots/figure5E.svg", width = 4, height = 5)
ggplot(sites_count_m %>% group_by(IntegrationStatus) %>% 
         mutate(Percentage = round(value / sum(value), 2)), 
       aes(x = IntegrationStatus, y = Percentage, fill = Region, cumulative = TRUE)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = function(x) paste0(x * 100)) +
  geom_col() +
  geom_text(aes(label = paste0(value, " (", Percentage * 100, "%)")), position = position_stack(vjust = 0.5)) +
  labs(title = "", x = "", y = "Count (Percentage)", fill = " ") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.background = element_rect(fill = NA),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", linewidth = 0.8, fill = NA))
dev.off()

# Gene expression z-score analysis
geneExpr_df <- read.csv("./data/all_genes_z_scores_no_batch_removal.csv", check.names = FALSE)
rownames(geneExpr_df) <- geneExpr_df$gene
geneExpr_df$gene <- NULL
z_df_102 <- geneExpr_df[, sample_102]
z_df_long <- z_df_102 %>% 
  tibble::rownames_to_column(var = "ConfidentGene.x") %>%
  pivot_longer(cols = -ConfidentGene.x, names_to = "SampleID", values_to = "value")
classified_wgs_sites <- classified_wgs_sites %>%
  left_join(z_df_long, by = c("SampleID", "ConfidentGene.x")) %>%
  rename(z_score = value)
classified_wgs_sites$abs_z <- abs(classified_wgs_sites$z_score)

# Violin plot for z-score (Figure 5F left)
svg("./plots/figure5F_left.svg", width = 3, height = 5)
ggplot(classified_wgs_sites, aes(x = group, y = z_score)) +
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("darkorange1", "skyblue")) +
  labs(title = " ", x = "") +
  scale_y_continuous(limits = c(-3, 10), breaks = seq(-3, 10, by = 3)) +
  geom_jitter(shape = 16, position = position_jitter(0.08), size = 1) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2.5, fill = "white") +
  stat_compare_means(comparisons = list(c("productive", "silent")), method = "t.test") +
  labs(y = "z-score of annotated gene") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.1, linetype = "dashed"),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.1, linetype = "dashed"),
        panel.border = element_rect(linetype = "solid", size = 0.8, fill = NA))
dev.off()

# Violin plot for absolute z-score (Figure 5F right)
svg("./plots/figure5F_right.svg", width = 3, height = 5)
ggplot(classified_wgs_sites, aes(x = group, y = abs_z)) +
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("darkorange1", "skyblue")) +
  labs(title = " ", x = "") +
  scale_y_continuous(limits = c(-3, 10), breaks = seq(-3, 10, by = 3), position = "right") +
  geom_jitter(shape = 16, position = position_jitter(0.08), size = 1) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2.5, fill = "white") +
  stat_compare_means(comparisons = list(c("productive", "silent")), method = "t.test") +
  labs(y = "|z-score| of annotated gene") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.1, linetype = "dashed"),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.1, linetype = "dashed"),
        panel.border = element_rect(linetype = "solid", size = 0.8, fill = NA))
dev.off()
