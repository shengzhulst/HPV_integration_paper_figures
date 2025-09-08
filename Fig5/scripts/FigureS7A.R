# Author: Shaomiao Xia
# Date: [Add date here]
# Description: Analysis of HPV integration status and immunonco scores, with violin plot visualization.

# Load required libraries
library(ComplexHeatmap)

# Set working directory
setwd("/Users/shaomiao/Desktop/Umich/SartorLab/HPV_integration/code_for_submission")

# Read sample status and sample lists
status <- read.csv("./data/all_HPV_samples_int_status_two_judgement_manullay.csv")
all_samples <- read.delim("./data/samples_105.txt", header=FALSE)$V1
DNA_pos_samples <- read.delim("./data/DNA_pos.txt", header=FALSE)$V1

# Remove unwanted samples and standardize sample names
remove_samples <- c("GS18055T", "GS18077T", "GS18107T")
all_samples <- gsub("-", ".", all_samples[!all_samples %in% remove_samples])
DNA_pos_samples <- gsub("-", ".", DNA_pos_samples[!DNA_pos_samples %in% remove_samples])
status$sample_name <- gsub("-", ".", status$sample_name)
status$sample_name <- ifelse(grepl("^\\d", status$sample_name), paste0("X", status$sample_name), status$sample_name)

# Identify DNA negative samples
DNA_neg_samples <- setdiff(all_samples, DNA_pos_samples)

# Subset status for relevant samples
samples_subset_status <- status[status$sample_name %in% all_samples,]

# Group samples by HPV integration status
status_groups <- list(
  R_pos_conf = samples_subset_status$sample_name[samples_subset_status$HPV_int_status_final == "HPVint(+) confident"],
  R_pos_like = samples_subset_status$sample_name[samples_subset_status$HPV_int_status_final == "HPVint(+) likely"],
  R_neg_conf = samples_subset_status$sample_name[samples_subset_status$HPV_int_status_final == "HPVint(-) confident"],
  R_neg_like = samples_subset_status$sample_name[samples_subset_status$HPV_int_status_final == "HPVint(-) likely"]
)

# Define sample groups
DpRp_samples_new <- intersect(DNA_pos_samples, c(status_groups$R_pos_conf, status_groups$R_pos_like))
DpRn_samples_new <- intersect(DNA_pos_samples, c(status_groups$R_neg_conf, status_groups$R_neg_like))
DnRp_samples_new <- intersect(DNA_neg_samples, c(status_groups$R_pos_conf, status_groups$R_pos_like))
DnRn_samples_new <- intersect(DNA_neg_samples, c(status_groups$R_neg_conf, status_groups$R_neg_like))

### Figure S7A
# Load z-score matrix and DEGs
z_df <- read.csv("./data/all_genes_z_scores_no_batch_removal_unfiltered_epsilon_1e-6.csv")
rownames(z_df) <- z_df$X
z_df$X <- NULL
qlf_ranked_sig <- read.csv("./data/qlf_ranked_sig_DEGs.csv")
top_DEGs <- intersect(qlf_ranked_sig$genes, rownames(z_df))

# Subset z-score matrix for top DEGs and relevant samples
z_df_all_samples <- z_df[top_DEGs, intersect(colnames(z_df), all_samples)]
DpRp_mat <- z_df[top_DEGs, DpRp_samples_new]
DpRn_mat <- z_df[top_DEGs, DpRn_samples_new]
DnRn_mat <- z_df[top_DEGs, DnRn_samples_new]

# Combine matrices and create annotation vectors
split_vec <- c(rep("DNA+RNA+", length(colnames(DpRp_mat))),
               rep("DNA+RNA-", length(colnames(DpRn_mat))),
               rep("DNA-RNA-", length(colnames(DnRn_mat))))
all_mat <- cbind(DpRp_mat, DpRn_mat, DnRn_mat)
sample_ids <- colnames(all_mat)
cohort_anno <- ifelse(startsWith(sample_ids, "GS"), "HVC", "TCGA")

# Load and process metadata
samples_meta <- read.csv("./data/all_scores_meta.csv")
samples_meta$Mixture <- gsub("-", ".", samples_meta$Mixture)
samples_meta <- samples_meta %>%
  mutate(Mixture = ifelse(startsWith(Mixture, "GS"), paste0(Mixture, "T"), Mixture)) %>%
  mutate(Mixture = ifelse(startsWith(Mixture, "TCGA"), paste0(Mixture, ".01A"), Mixture))
samples_meta$Mixture <- gsub("TCGA.BA.4077.01A", "TCGA.BA.4077.01B", samples_meta$Mixture)
rownames(samples_meta) <- samples_meta$Mixture
samples_meta$Mixture <- NULL
meta_sub_df <- samples_meta[colnames(all_mat),]
meta_sub_df <- meta_sub_df %>%
  mutate(Sex = ifelse(Sex_F == "FALSE", "Male", "Female")) %>%
  mutate(Subtype = ifelse(Nosubtype_IMU == "FALSE", "KRT", "IMU"))

# Prepare ImmunOnco scores for annotation
rownames(immunonco) <- immunonco$Mixture
escores_sub_df <- immunonco[colnames(all_mat),]

# Plot heatmap and save as SVG
svg("./plots/figureS7a.svg", width = 7, height = 5)
Heatmap(all_mat, name = "z_scores",
        column_title = "Samples",
        row_title = "Differentially Expressed Genes",
        column_labels = rep("", length(split_vec)),
        row_labels = rep("", nrow(all_mat)),
        column_title_side = "bottom",
        top_annotation = columnAnnotation(
          Category = split_vec,
          Cohorts = cohort_anno,
          Sex = meta_sub_df$Sex,
          ImmunOncoScores = anno_points(escores_sub_df$e1.e2.e5.l1.l2__avg, size = unit(1.5, "mm"), extend = 0.1),
          col = list(
            Category = c("DNA+RNA+" = "darkorange1", "DNA+RNA-" = "skyblue", "DNA-RNA-" = "cornflowerblue"),
            Cohorts = c("HVC" = "grey", "TCGA" = "black"),
            Sex = c("Male" = "lightgreen", "Female" = "pink")
          ),
          simple_anno_size = unit(0.35, "cm")
        ))
dev.off()
