########README
#####integrate the ATAC and RNA, then using the statistical 
#ranking (not better than RF, Xgboost and RL with ranking gene IP6K1 top )

### in another one, will use ML for ranking

#-----------------------------------------------
# STEP 0: Setup and Installation
#-----------------------------------------------

# Define required packages for the complete pipeline
required_bioc <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", 
                   "clusterProfiler", "ReactomePA", "enrichplot", "AnnotationDbi",
                   "GenomicRanges", "rtracklayer", "EnrichedHeatmap", 
                   "ComplexHeatmap", "DESeq2")

required_cran <- c("data.table", "ggplot2", "pheatmap", "RColorBrewer", 
                   "ggrepel", "ggpubr", "msigdbr", "randomForest", "xgboost")

# 1. Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Using update=FALSE and ask=FALSE to prevent session hangs in 2025
new_bioc <- required_bioc[!(required_bioc %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# 2. Install missing CRAN packages
new_cran <- required_cran[!(required_cran %in% installed.packages()[,"Package"])]
if(length(new_cran)) install.packages(new_cran)

# 3. Load all core libraries required for integration and enrichment
library(data.table)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(ReactomePA)   # Specifically added for pathway analysis
library(enrichplot)   # Required for dotplot()
library(msigdbr)      # 2025 standard for accessing pathway databases

# Set global TxDb for the hg38 genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#-----------------------------------------------
# STEP 1: Load and Filter Data
#-----------------------------------------------
# 1.1: RNA-seq
rna_data <- fread("C:/Users/qiqi5/Downloads/m01/lei-m/ATACseq_RNAseq_diff/RNA_seq_differential_expression.tsv")
sig_genes <- rna_data[abs(log2FoldChange) > 1 & padj < 0.05]

# 1.2: ATAC-seq
atac_diff <- fread("C:/Users/qiqi5/Downloads/m01/lei-m/ATACseq_RNAseq_diff/ATAC_seq_differential_peaks.tsv")
sig_atac <- atac_diff[abs(log2FoldChange) > 1 & padj < 0.05]

# Convert to GRanges
sig_atac_gr <- makeGRangesFromDataFrame(sig_atac, keep.extra.columns = TRUE,
                                        seqnames.field = "chr", start.field = "start", end.field = "end")

#-----------------------------------------------
# STEP 2: Annotation
#-----------------------------------------------
# Annotate peaks
diff_peak_anno <- annotatePeak(sig_atac_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Hs.eg.db")

# Convert to data.table
diff_peak_anno_df <- as.data.table(as.data.frame(diff_peak_anno))

# Add Gene Symbols using AnnotationDbi::mapIds
diff_peak_anno_df[, gene_symbol := mapIds(org.Hs.eg.db, keys = as.character(geneId), 
                                          column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")]

# Visualization
pdf("ATAC_diffpeak_plots.pdf", width = 10, height = 8)
plotAnnoPie(diff_peak_anno)
plotDistToTSS(diff_peak_anno)
dev.off()

##################################################################################
##################################################################################
#-----------------------------------------------
# STEP 3: Integration
#-----------------------------------------------
# 3.1 Prepare for merge - Creating unique keys for peak matching
sig_atac[, peak_key := paste(chr, start, end, sep = "_")]
diff_peak_anno_df[, peak_key := paste(seqnames, start, end, sep = "_")]

# Select necessary columns and merge accessibility info
atac_with_genes <- diff_peak_anno_df[, .(peak_key, annotation, gene_symbol, geneId, distanceToTSS)]
atac_with_genes <- merge(atac_with_genes, sig_atac[, .(peak_key, log2FoldChange, padj)], by = "peak_key")
setnames(atac_with_genes, old = c("log2FoldChange", "padj"), new = c("atac_log2FC", "atac_padj"))

# Prepare RNA data
rna_sig <- sig_genes[, .(gene_name, log2FoldChange, padj)]
setnames(rna_sig, old = c("log2FoldChange", "padj"), new = c("rna_log2FC", "rna_padj"))

# 3.2 Final Integration (Inner Join)
integrated_data <- merge(atac_with_genes, rna_sig, by.x = "gene_symbol", by.y = "gene_name", allow.cartesian = TRUE)

# 3.3 Categorization and Statistics
integrated_data[, concordance := fcase(
  atac_log2FC > 0 & rna_log2FC > 0, "Concordant up",
  atac_log2FC < 0 & rna_log2FC < 0, "Concordant down",
  atac_log2FC > 0 & rna_log2FC < 0, "Discordant (open/down)",
  atac_log2FC < 0 & rna_log2FC > 0, "Discordant (closed/up)",
  default = "Other"
)]

# Output Summary
correlation <- cor.test(integrated_data$atac_log2FC, integrated_data$rna_log2FC)
cat("\n--- Integration Results ---\n")
cat("Integrated pairs:", nrow(integrated_data), "\n")
cat("Correlation:", round(correlation$estimate, 3), "(p =", format.pval(correlation$p.value), ")\n")
print(integrated_data[, .N, by = concordance])


#####################################################
#####################################################
#-----------------------------------------------
# STEP 4: Visualize the Integrated Results
#-----------------------------------------------

#=============================================
# 4.1: Calculate Integrated Scores for Ranking
#=============================================
# We use fcoalesce to handle potential NAs and pmin to cap p-values for scoring
integrated_data[, `:=` (
  # Standard product of fold changes
  combined_score = atac_log2FC * rna_log2FC,
  
  # Z-score normalized product (gives equal weight to both data types)
  z_combined_score = scale(atac_log2FC)[,1] * scale(rna_log2FC)[,1],
  
  # Statistical significance score (using a small offset to avoid log10(0))
  sig_score = -log10(pmax(atac_padj, 1e-300)) - log10(pmax(rna_padj, 1e-300)),
  
  # Weighted score: Magnitude * Significance
  weighted_score = (atac_log2FC * rna_log2FC) * (-log10(pmax(atac_padj * rna_padj, 1e-300))),
  
  # Distance-weighted score (exponential decay: weight drops as distance increases)
  distance_weight = exp(-pmin(abs(distanceToTSS), 100000)/10000)
)]

integrated_data[, distance_weighted_score := combined_score * distance_weight]

#=============================================
# 4.2: Integrated Scatterplot (ATAC vs RNA)
#=============================================
scatter_plot <- ggplot(integrated_data, 
                       aes(x = atac_log2FC, y = rna_log2FC, 
                           color = concordance,
                           size = -log10(pmax(atac_padj * rna_padj, 1e-300)))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("Concordant up" = "#E63946", 
                                "Concordant down" = "#1D3557",
                                "Discordant (open/down)" = "#F4A261", 
                                "Discordant (closed/up)" = "#7B2CBF",
                                "Other" = "gray70")) +
  scale_size_continuous(range = c(1, 6), name = "-log10(Combined P)") +
  labs(x = "ATAC-seq log2FC", y = "RNA-seq log2FC",
       title = "Chromatin Accessibility vs. Gene Expression",
       subtitle = paste("Pearson r:", round(correlation$estimate, 3), 
                        "| p:", format.pval(correlation$p.value))) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  # Label top 15 most "impactful" genes based on weighted score
  geom_text_repel(data = integrated_data[order(-abs(weighted_score))][1:15],
                  aes(label = gene_symbol), size = 3, color = "black",
                  max.overlaps = 30, box.padding = 0.5)

ggsave("ATAC_RNA_scatter.pdf", scatter_plot, width = 8, height = 6)

#=============================================
# 4.3: Integration Volcano Plot
#=============================================
integrated_data[, neg_log10_pval := -log10(pmax(atac_padj * rna_padj, 1e-300))]

volcano_plot <- ggplot(integrated_data, 
                       aes(x = combined_score, y = neg_log10_pval, color = concordance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Concordant up" = "#E63946", 
                                "Concordant down" = "#1D3557",
                                "Discordant (open/down)" = "#F4A261", 
                                "Discordant (closed/up)" = "#7B2CBF",
                                "Other" = "gray70")) +
  theme_minimal() +
  labs(x = "Combined Score (ATAC log2FC * RNA log2FC)", 
       y = "-log10(Combined P-value)",
       title = "Integrated Significance Volcano") +
  geom_text_repel(data = integrated_data[order(-neg_log10_pval)][1:15],
                  aes(label = gene_symbol), size = 3, color = "black")

ggsave("Integration_volcano.pdf", volcano_plot, width = 8, height = 6)

#=============================================
# 4.4: Peak-Gene Distance Analysis
#=============================================
integrated_data[, distance_kb := abs(distanceToTSS) / 1000]

distance_hist <- ggplot(integrated_data, aes(x = distance_kb, fill = concordance)) +
  geom_histogram(bins = 40, alpha = 0.8, color = "white", linewidth = 0.1) +
  facet_wrap(~ concordance, scales = "free_y") +
  scale_fill_manual(values = c("Concordant up" = "#E63946", 
                               "Concordant down" = "#1D3557",
                               "Discordant (open/down)" = "#F4A261", 
                               "Discordant (closed/up)" = "#7B2CBF",
                               "Other" = "gray70")) +
  theme_light() +
  labs(x = "Distance from TSS (kb)", y = "Peak Count",
       title = "Genomic Distribution by Concordance") +
  theme(legend.position = "none")

ggsave("Distance_histograms.pdf", distance_hist, width = 10, height = 7)


#-----------------------------------------------
# STEP 4.5: Statistical Ranking (Non-ML)
#-----------------------------------------------

# 1. Calculate a Combined Evidence Score
# Logic: (ATAC Strength * RNA Strength) / (Distance Penalty)
# We use abs() because both 'Concordant Up' and 'Concordant Down' are important.

integrated_data[, stat_rank_score := (abs(atac_log2FC) * abs(rna_log2FC)) * 
                  (-log10(pmax(atac_padj * rna_padj, 1e-300))) * 
                  distance_weight]

# 2. Handle multiple peaks per gene (Collapse to the best regulatory event)
# This gives us one row per gene, keeping the highest score found.
final_stat_ranking <- integrated_data[, .(
  max_stat_score = max(stat_rank_score),
  concordance = concordance[which.max(stat_rank_score)],
  atac_fc = atac_log2FC[which.max(stat_rank_score)],
  rna_fc = rna_log2FC[which.max(stat_rank_score)],
  dist_to_tss = distanceToTSS[which.max(stat_rank_score)]
), by = gene_symbol][order(-max_stat_score)]

# 3. View the Top 10 Statistically Ranked Genes
print(head(final_stat_ranking, 10))

# 4. Export for GO Enrichment
fwrite(final_stat_ranking, "Statistical_Ranked_Genes_NoML.csv")

######################################################################
#################### Survial part###############################

#-----------------------------------------------
# STEP 9: Advanced Survival Analysis (KM, Cox, RSF, XGBoost)
#-----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install both at once
BiocManager::install(c("randomForestSRC", "xgboost"))

library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(survival)
library(survminer)
library(randomForestSRC)
library(xgboost)

###############################################################
######################################################################
#################### Survial part###############################
library(UCSCXenaTools)
library(dplyr)

# 1. Query the full IlluminaHiSeq RNA-seq dataset for LUAD
luad_hiseq_info <- XenaData %>% 
  dplyr::filter(XenaCohorts == "TCGA Lung Adenocarcinoma (LUAD)", 
                XenaDatasets == "TCGA.LUAD.sampleMap/HiSeqV2") # Standard RNA-seq matrix

# 2. Download the full expression matrix
luad_hiseq <- luad_hiseq_info %>%
  XenaGenerate() %>%
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()


library(data.table)

# Convert to data.table if not already
hiseq_dt <- as.data.table(luad_hiseq)

# 1. Find the column that contains gene names
# Usually, it's the first column. Let's check for IP6K1 there.
gene_col <- colnames(hiseq_dt)[1]

# Try to find IP6K1 or Entrez ID 9807
target_row <- hiseq_dt[get(gene_col) == "IP6K1" | get(gene_col) == "9807"]



################################################################




# Prepare data objects
surv_df <- as.data.table(luad_data$LUAD_survival.txt)
hiseq_dt <- as.data.table(luad_data$HiSeqV2)
clinical_dt <- as.data.table(luad_data$LUAD_clinicalMatrix)

target_gene <- "IP6K1"
gene_col <- colnames(hiseq_dt)[1]
target_row <- hiseq_dt[get(gene_col) == target_gene]

if (nrow(target_row) > 0) {
  # 1. Format Expression
  gene_expr <- t(target_row[, -1, with = FALSE])
  expr_clean <- data.table(sample = rownames(gene_expr), 
                           Expression = as.numeric(gene_expr[,1]))
  
  # 2. Merge All Layers
  meta_data <- clinical_dt[, .(sampleID, 
                               age = age_at_initial_pathologic_diagnosis, 
                               stage = pathologic_stage)]
  
  plot_data <- merge(surv_df, expr_clean, by = "sample")
  plot_data <- merge(plot_data, meta_data, by.x = "sample", by.y = "sampleID")
  
  # 3. CRITICAL DATA CLEANING (Fixes Factor Error)
  # Simplify Stage and convert to Factor
  plot_data[, stage_simple := sub(" [AB-C]$| [AB-C][1-2]$| [1-2]$|A$|B$", "", stage)]
  plot_data <- plot_data[stage_simple %in% c("Stage I", "Stage II", "Stage III", "Stage IV")]
  plot_data[, stage_simple := as.factor(stage_simple)]
  
  # Ensure Age is numeric and remove NAs
  plot_data[, age := as.numeric(age)]
  plot_data <- na.omit(plot_data, cols = c("Expression", "OS.time", "age", "stage_simple"))
  
  setnames(plot_data, "Expression", target_gene)
  
  # --- 9.3: Kaplan-Meier ---
  median_val <- median(plot_data[[target_gene]])
  plot_data[, Group := ifelse(get(target_gene) > median_val, "High", "Low")]
  fit_km <- survfit(Surv(OS.time, OS) ~ Group, data = plot_data)
  print(ggsurvplot(fit_km, data = plot_data, pval = TRUE, risk.table = TRUE,
                   title = paste("KM Plot:", target_gene)))
  
  # --- 9.4: Multivariable Cox PH Model ---
  # Using explicit column names to avoid 'get' environment issues
  cox_model <- coxph(Surv(OS.time, OS) ~ IP6K1 + age + stage_simple, data = plot_data)
  message("\n--- Multivariable Cox Model Results ---")
  print(summary(cox_model))
  print(cox.zph(cox_model))
  
  # --- 9.5: Random Survival Forest (RSF) ---
  set.seed(2025)
  rsf_model <- rfsrc(Surv(OS.time, OS) ~ IP6K1 + age + stage_simple, 
                     data = plot_data, 
                     ntree = 500)
  message("\n--- RSF Variable Importance (VIMP) ---")
  print(rsf_model$importance)
  
  # --- 9.6: XGBoost Survival ---
  # Convert factors to dummy variables
  xgb_features <- model.matrix(Surv(OS.time, OS) ~ IP6K1 + age + stage_simple, data = plot_data)[,-1]
  xgb_label <- ifelse(plot_data$OS == 1, plot_data$OS.time, -plot_data$OS.time)
  
  xgb_surv <- xgboost(data = xgb_features, label = xgb_label, 
                      params = list(objective = "survival:cox", learning_rate = 0.05),
                      nrounds = 50, verbosity = 0)
  
  importance_surv <- xgb.importance(model = xgb_surv)
  message("\n--- XGBoost Survival Importance ---")
  print(importance_surv)
  
} else {
  message("Gene not found.")
}



















# --- 9.1: Data Integration (Clinical + Survival + Expression) ---
surv_df <- as.data.table(luad_data$LUAD_survival.txt)
hiseq_dt <- as.data.table(luad_data$HiSeqV2)
clinical_dt <- as.data.table(luad_data$LUAD_clinicalMatrix)

# Extract IP6K1
target_gene <- "IP6K1"
gene_col <- colnames(hiseq_dt)[1]
target_row <- hiseq_dt[get(gene_col) == target_gene]

if (nrow(target_row) > 0) {
  # Format Expression
  gene_expr <- t(target_row[, -1, with = FALSE])
  expr_clean <- data.table(sample = rownames(gene_expr), 
                           Expression = as.numeric(gene_expr[,1]))
  
  # Merge All Layers
  # We select Age and Pathologic Stage from the clinical matrix
  meta_data <- clinical_dt[, .(sampleID, 
                               age = age_at_initial_pathologic_diagnosis, 
                               stage = pathologic_stage)]
  

  plot_data <- merge(surv_df, expr_clean, by = "sample")
  plot_data <- merge(plot_data, meta_data, by.x = "sample", by.y = "sampleID")
  
  
  # Rename Expression to Gene Name for model clarity
  setnames(plot_data, "Expression", target_gene)
  
  # --- 9.3: Kaplan-Meier (Visual) ---
  median_val <- median(plot_data[[target_gene]])
  plot_data[, Group := ifelse(get(target_gene) > median_val, "High", "Low")]
  fit_km <- survfit(Surv(OS.time, OS) ~ Group, data = plot_data)
  
  print(ggsurvplot(fit_km, data = plot_data, pval = TRUE, risk.table = TRUE,
                   title = paste("KM Plot:", target_gene)))
  
  # --- 9.4: Multivariable Cox PH Model ---
  # Standardized for 2025: Does the gene work independent of Age and Stage?
  cox_model <- coxph(Surv(OS.time, OS) ~ get(target_gene) + age + stage_simple, data = plot_data)
  message("\n--- Multivariable Cox Model Results ---")
  print(summary(cox_model))
  
  # Check PH Assumption
  print(cox.zph(cox_model))
  
  # --- 9.5: Random Survival Forest (ML) ---
  # RSF handles non-linear interactions between Age/Stage/Gene
  set.seed(2025)
  rsf_model <- rfsrc(Surv(OS.time, OS) ~ ., 
                     data = plot_data[, .(OS.time, OS, IP6K1, age, stage_simple)], 
                     ntree = 500)
  
  message("\n--- RSF Variable Importance (VIMP) ---")
  print(rsf_model$importance)
  plot(rsf_model)
  
  # --- 9.6: XGBoost Survival (Gradient Boosting) ---
  # Prep matrix
  xgb_features <- model.matrix(~ IP6K1 + age + stage_simple - 1, data = plot_data)
  # Target for Cox objective: >0 for events, <0 for censored
  xgb_label <- ifelse(plot_data$OS == 1, plot_data$OS.time, -plot_data$OS.time)
  
  xgb_surv <- xgboost(data = xgb_features, label = xgb_label, 
                      params = list(objective = "survival:cox", learning_rate = 0.05),
                      nrounds = 50, verbosity = 0)
  
  # Feature Importance
  importance_surv <- xgb.importance(model = xgb_surv)
  message("\n--- XGBoost Survival Importance ---")
  print(importance_surv)
  
} else {
  message("Gene not found.")
}
