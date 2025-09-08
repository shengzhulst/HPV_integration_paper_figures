# Author: Shaomiao Xia
# Description: Analysis of HPV integration status and immunonco scores, with violin plot visualization Figure 5G.

# Load required libraries
library(ggpubr)
library(vcd)

# Set working directory
setwd("HPV_integration_paper_figures/Fig5")

# Read sample status and sample lists
status <- read.csv("./data/all_HPV_samples_int_status_two_judgement_manullay.csv")
all_samples <- read.delim("./data/samples_105.txt", header = FALSE)$V1
DNA_pos_samples <- read.delim("./data/DNA_pos.txt", header = FALSE)$V1

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

# Read immunonco scores and standardize mixture names
immunonco <- read.csv("./data/scRNA_all_int_cohort_adj_100Constant.csv")
immunonco$Mixture <- gsub("EGAF[0-9]+_", "", immunonco$Mixture)
immunonco$Mixture <- gsub("-", ".", immunonco$Mixture)
immunonco$Mixture <- ifelse(grepl("^\\d", immunonco$Mixture), paste0("X", immunonco$Mixture), immunonco$Mixture)

# Assign group labels for plotting
DpRp_e_score <- immunonco[immunonco$Mixture %in% DpRp_samples_new,]
DpRn_e_score <- immunonco[immunonco$Mixture %in% DpRn_samples_new,]
DnRn_e_score <- immunonco[immunonco$Mixture %in% DnRn_samples_new,]

DpRp_e_score$group <- "DNA+RNA+"
DpRn_e_score$group <- "DNA+RNA-"
DnRn_e_score$group <- "DNA-RNA-"

# Prepare data for violin plot
# ADJ_e1.e2.e5.l1.l2__avg: Adjusted ImmunOnco scores
# group: Sample group
df_violin_immuno <- rbind(
  DpRp_e_score[c("ADJ_e1.e2.e5.l1.l2__avg", "group")],
  DpRn_e_score[c("ADJ_e1.e2.e5.l1.l2__avg", "group")],
  DnRn_e_score[c("ADJ_e1.e2.e5.l1.l2__avg", "group")]
)

# Violin plot of HPV Gene Ratio by sample group
svg("./plots/figure5G.svg", width = 5, height = 5)
p <- ggplot(df_violin_immuno, aes(x = group, y = ADJ_e1.e2.e5.l1.l2__avg)) +
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("cornflowerblue", "skyblue", "darkorange1")) +
  labs(title = " ", x = "", y = "HPV Gene Ratio") +
  scale_y_continuous(limits = c(-4, 6), breaks = c(-4, -2, 0, 2, 4, 6)) +
  geom_jitter(shape = 16, position = position_jitter(0.05), size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white") +
  stat_compare_means(comparisons = list(c("DNA+RNA+", "DNA+RNA-"), c("DNA+RNA+", "DNA-RNA-"), c("DNA+RNA-", "DNA-RNA-")), test = "t.test") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.1, linetype = "dashed"),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.1, linetype = "dashed"),
        panel.border = element_rect(linetype = "solid", size = 0.8, fill = NA))
print(p)
dev.off()
