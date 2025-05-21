# Load required libraries
library(vcfR)
library(adegenet)
library(ggplot2)
library(viridis)
library(gridExtra)
library(reshape2)
library(dplyr)

# Step 1: Load VCF data
setwd("~/Affinis/affinis_FULL/affinis_strict")

# vcf <- read.vcfR("affinis_full_50.recode.vcf")  # 59,021 variants
# vcf <- read.vcfR("affinis_full_75.recode.vcf")  # 31,233 variants
#vcf <- read.vcfR("affinis_full_90.recode.vcf")  # 10,991 variants
vcf <- read.vcfR("affinis_full_80.recode.vcf")  # 24,863 variants

# Step 2: Convert to genlight
geno <- vcfR2genlight(vcf)

# Step 3: Load and clean metadata
metadata <- read.csv("affinis_metadata_floricomplex.csv", header = TRUE, sep = ",")
head(metadata)
metadata$SampleID <- as.character(metadata$Sample_ID)

# Step 4: Match and assign population info to genlight
gl_samples <- indNames(geno)
metadata <- metadata[metadata$SampleID %in% gl_samples, ]
metadata <- metadata[match(gl_samples, metadata$SampleID), ]  # Ensure same order

# Double check no NA assignments
stopifnot(all(metadata$SampleID == gl_samples))  # Safety check
stopifnot(!any(is.na(metadata$Taxon)))           # No NAs allowed in Taxon

# Assign populations to genlight object
pop(geno) <- factor(metadata$Taxon)

# Step 5: Find optimal number of PCs (if needed)
set.seed(123456)
# dapc_res <- dapc(geno, n.pca = 15, n.da = 4)
# opt_pca <- optim.a.score(dapc_res, n.pca = 1:15, n.da = 4, n.sim = 1000)
# print(opt_pca)  # <- visually choose best n.pca
### 7 PC axes for affinis90strict

# Step 6: DAPC + posterior plotting function
library(adegenet)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(gridExtra)
library(purrr)
library(combinat)

perform_multi_dapc <- function(geno, metadata,
                               pca_values = c(2, 4, 5, 7),
                               n_da = 2,
                               label_prefix = "Group",
                               exclude_taxa = NULL,
                               save_plots = TRUE,
                               output_dir = "dapc_results") {
  dir.create(output_dir, showWarnings = FALSE)
  
  # Filter metadata and exclude specified taxa
  metadata <- metadata %>%
    filter(Sample_ID %in% indNames(geno)) %>%
    filter(!(Taxon %in% exclude_taxa)) %>%
    group_by(Taxon) %>% filter(n() > 1) %>% ungroup() %>%
    filter(complete.cases(.))
  rownames(metadata) <- metadata$Sample_ID
  
  # Match geno to metadata
  geno <- geno[indNames(geno) %in% metadata$Sample_ID]
  geno <- geno[match(metadata$Sample_ID, indNames(geno))]
  
  # Set population grouping
  pop(geno) <- as.factor(metadata$Taxon)
  
  results_list <- list()
  
  for (pca_num in pca_values) {
    cat("Running DAPC with", pca_num, "PCs...\n")
    
    tryCatch({
      dapc_result <- dapc(geno, n.pca = pca_num, n.da = n_da)
      
      df_scores <- as.data.frame(dapc_result$ind.coord)
      df_scores$Sample_ID <- rownames(df_scores)
      df_scores$Taxon <- metadata$Taxon
      
      # Determine available LD axes
      ld_axes <- colnames(df_scores)[grepl("^LD", colnames(df_scores))]
      ld_combos <- combn(ld_axes, 2, simplify = FALSE)
      
      # Create scatter plots for each LD pair
      scatter_plots <- list()
      for (combo in ld_combos) {
        ld_x <- combo[1]
        ld_y <- combo[2]
        
        # Convex hulls for plotting
        hull_data <- df_scores %>%
          group_by(Taxon) %>%
          slice(chull(.data[[ld_x]], .data[[ld_y]]))
        
        scatter_plot <- ggplot(df_scores, aes_string(x = ld_x, y = ld_y, color = "Taxon")) +
          geom_polygon(data = hull_data, aes_string(fill = "Taxon", group = "Taxon"),
                       alpha = 0.15, color = NA) +
          geom_point(size = 3) +
          geom_text_repel(aes(label = Sample_ID), size = 3, max.overlaps = 10) +
          theme_minimal() +
          labs(title = paste0("DAPC: ", ld_x, " vs ", ld_y, " (PCs=", pca_num, ")"),
               x = ld_x, y = ld_y) +
          theme(legend.position = "right")
        
        scatter_plots[[paste0(ld_x, "_vs_", ld_y)]] <- scatter_plot
        
        # Save plot
        if (save_plots) {
          filename <- file.path(output_dir,
                                paste0("scatter_", ld_x, "_vs_", ld_y, "_PCs", pca_num, ".png"))
          ggsave(filename, plot = scatter_plot, width = 8, height = 6)
        }
      }
      
      # Posterior probabilities
      posterior_df <- as.data.frame(dapc_result$posterior)
      posterior_df$Sample_ID <- rownames(posterior_df)
      posterior_df$Taxon <- metadata$Taxon
      
      posterior_long <- posterior_df %>%
        pivot_longer(-c(Sample_ID, Taxon), names_to = "Group", values_to = "Posterior")
      
      posterior_plot <- ggplot(posterior_long, aes(x = Sample_ID, y = Posterior, fill = Group)) +
        geom_col(position = "stack") +
        facet_wrap(~Taxon, scales = "free_x", nrow = 1) +
        theme_minimal() +
        labs(title = paste0("Posterior Probabilities - ", label_prefix, " (PCs=", pca_num, ")"),
             x = "Sample ID", y = "Posterior Probability") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
      
      if (save_plots) {
        ggsave(file.path(output_dir, paste0("posterior_PCs", pca_num, ".png")),
               plot = posterior_plot, width = 10, height = 5)
      }
      
      results_list[[paste0("PC", pca_num)]] <- list(
        dapc = dapc_result,
        scatter_plots = scatter_plots,
        posterior_plot = posterior_plot
      )
      
    }, error = function(e) {
      message(sprintf("⚠️ Skipped PCA = %d due to error: %s", pca_num, e$message))
    })
  }
  
  # Composite saving
  if (save_plots && length(results_list) > 0) {
    for (pc_key in names(results_list)) {
      scatters <- results_list[[pc_key]]$scatter_plots
      if (length(scatters) > 0) {
        ggsave(file.path(output_dir, paste0("composite_", pc_key, ".png")),
               plot = marrangeGrob(grobs = scatters, nrow = 2, ncol = 2), width = 12, height = 8)
      }
      ggsave(file.path(output_dir, paste0("posterior_", pc_key, ".png")),
             plot = results_list[[pc_key]]$posterior_plot, width = 10, height = 5)
    }
  }
  
  return(results_list)
}


# Step 7: Run DAPC + plot
results <- perform_multi_dapc(
  geno = geno,
  metadata = metadata,
  pca_values = c(2, 4, 5, 7),
  n_da = 4,
  label_prefix = "Affinis_strict90",
  #exclude_taxa = c("Rosacea", "Impostor"),
  save_plots = TRUE,
  output_dir = "Affinis_DAPC_Output_strict90"
)

results
# Step 8: View results
grid.arrange(grobs = results$scatter_plots[1:4], ncol = 2)
print(results$posterior_plot)

# Step 9: Save plots

#Optional: create specific output directory
output_dir <- "dapc_plots_50"
dir.create(output_dir, showWarnings = FALSE)

## Save each scatter plot individually
for (i in seq_along(results$scatter_plots)) {
  plot_name <- paste0(output_dir, "/scatter_LD_plot_", i, ".png")
  ggsave(filename = plot_name, plot = results$scatter_plots[[i]],
         width = 6, height = 5, dpi = 300)
}

## Save composite
library(gridExtra)

composite_plot <- grid.arrange(grobs = results$scatter_plots, ncol = 2)

ggsave(filename = file.path(output_dir, "scatter_LD_composite.png"),
       plot = composite_plot, width = 12, height = 10, dpi = 300)

## Save Posterior Probability Plot
# Save selected plots
ggsave("Viola_Affinis_DAPC_PC7_LD1_LD2.pdf", results$dapc_plots[["Affinis_strict50_nPCA7_LD1_LD2"]], width = 6, height = 5)
ggsave("Viola_Affinis_DAPC_PC7_LD1_LD3.pdf", results$dapc_plots[["Affinis_strict50_nPCA7_LD1_LD3"]], width = 6, height = 5)
ggsave("Viola_Affinis_Structure_PC7.pdf", results$structure_plots[["Affinis_strict50_nPCA7"]], width = 7, height = 3)
