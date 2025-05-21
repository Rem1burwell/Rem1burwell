###############################################################################
#################### Full affinis data ###########################################################
#################### Remington N. Burwell Ohio University #########################################################
library(LEA)        # For sNMF
library(ggplot2)    # For plotting
library(reshape2)   # To reshape the Q-matrix into long format for ggplot
library(ggrepel)    # Prevents overlapping text labels
library(gridExtra)  # For arranging multiple plots in a grid
library(tidyr)      # For data manipulation
library(dplyr)      # For data manipulation
library(vcfR)       # Handling vcf objects in R
library(adegenet)
library(ade4)
library(ape)
library(poppr)

# Step 1: Load VCF data
setwd("~/Desktop/Affinis_genomics")

vcf <- read.vcfR("affinis_full_90.recode.vcf")  # 10,991 variants
vcf2geno("affinis_full_90.recode.vcf") # generate geno file for sNMF
#project = load.snmfProject("affinis_full_90.recode.snmfProject")

#### Load in vcf and convert to geno file for sNMF
#vcf <- read.vcfR("affinis_full_80.recode.vcf") # 24863 snps
#vcf2geno("affinis_full_80.recode.vcf") # generate geno file for sNMF
### Load in sNMF project if started previously
#project = load.snmfProject("affinis_full_80.recode.snmfProject")

#### Run snmf with the LEA library testing K values 2-20 for 100 repetitions
obj <- snmf("affinis_full_90.recode.geno", K = 2:20, 
            project = "new", entropy = TRUE,
            repetitions = 100, ploidy = 10, CPU = 6, alpha = 10,
            tolerance = 0.00001, percentage = 0.1, 6200, seed = 6942069)

#### plot results of K
plot(obj, col = "blue", pch = 19)
###############################################################################
# ###############################################################################
#### Delta K (Evanno Method for determining K)
# Define K values to test
K_range <- 2:20 # Avoid K = 1 since ΔK is undefined

# Compute cross-entropy for each K
entropy_values <- sapply(K_range, function(k) min(cross.entropy(obj, K = k)))

# Compute ΔK (rate of change in cross-entropy)
delta_K <- abs(diff(entropy_values) / entropy_values[-length(entropy_values)])

# Convert to data frame
deltaK_df <- data.frame(K = K_range[-1], DeltaK = delta_K)

# Plot ΔK to find the best K
ggplot(deltaK_df, aes(x = K, y = DeltaK)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue", linetype = "dashed") +
  labs(title = "ΔK Plot for sNMF Model Selection",
       x = "Number of Clusters (K)",
       y = "ΔK (Rate of Change in Cross-Entropy)") +
  theme_minimal() +
  geom_vline(xintercept = deltaK_df$K[which.max(deltaK_df$DeltaK)], 
             color = "red", linetype = "dotted", linewidth = 1) +
  annotate("text", x = deltaK_df$K[which.max(deltaK_df$DeltaK)] + 0.5, 
           y = max(deltaK_df$DeltaK), label = "Best K", color = "red")
##### 2
###############################################################################
################################################################################
#### Cross Entropy Plots
# Compute cross-entropy values for each K
entropy_values <- sapply(K_range, function(k) min(cross.entropy(obj, K = k)))

# Convert to data frame for plotting
entropy_df <- data.frame(K = K_range, CrossEntropy = entropy_values)

# Find K with minimum cross-entropy
best_K_entropy <- entropy_df$K[which.min(entropy_df$CrossEntropy)]

# Plot Cross-Entropy vs. K
ggplot(entropy_df, aes(x = K, y = CrossEntropy)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue", linetype = "dashed") +
  labs(title = "Cross-Entropy for sNMF Models",
       x = "Number of Clusters (K)",
       y = "Cross-Entropy") +
  theme_minimal() +
  geom_vline(xintercept = best_K_entropy, 
             color = "red", linetype = "dotted", linewidth = 1) +
  annotate("text", x = best_K_entropy + 0.5, 
           y = min(entropy_df$CrossEntropy), label = "Best K", color = "red")
#### 10
###############################################################################
################################################################################
#### Mean Squared Error
# Extract cross-validation errors (using cross.entropy function)
MSE_values <- sapply(K_range, function(k) min(cross.entropy(obj, K = k)))
# Convert to data frame for ggplot
MSE_df <- data.frame(K = K_range, MSE = MSE_values)

# Find K with minimum MSE
best_K_MSE <- MSE_df$K[which.min(MSE_df$MSE)]

# Plot MSE vs. K
ggplot(MSE_df, aes(x = K, y = MSE)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue", linetype = "dashed") +
  labs(title = "Cross-Validation MSE for sNMF Models",
       x = "Number of Clusters (K)",
       y = "Mean Squared Error (MSE)") +
  theme_minimal() +
  geom_vline(xintercept = best_K_MSE, 
             color = "red", linetype = "dotted", linewidth = 1) +
  annotate("text", x = best_K_MSE + 0.5, 
           y = min(MSE_df$MSE), label = "Best K", color = "red")
### 10
##############################################################
#### K Value Table 
# Compute cross-entropy for each K
entropy_values <- sapply(K_range, function(k) min(cross.entropy(obj, K = k)))

# Compute ΔK (rate of change in cross-entropy)
delta_K <- abs(diff(entropy_values) / entropy_values[-length(entropy_values)])
deltaK_df <- data.frame(K = K_range[-1], DeltaK = delta_K)

# Compute MSE (same as cross-entropy in sNMF)
MSE_values <- entropy_values  

# Create a ranking table for K selection
ranking_table <- data.frame(
  K = K_range,
  CrossEntropy = entropy_values,
  MSE = MSE_values
)

# Add ΔK (set NA for K = 2 since it has no ΔK value)
ranking_table$DeltaK <- c(NA, delta_K)

# Rank K values for each method
ranking_table <- ranking_table %>%
  arrange(CrossEntropy) %>%  # Rank Cross-Entropy (low to high)
  mutate(Rank_CE = rank(CrossEntropy, ties.method = "first")) %>%
  arrange(MSE) %>%  # Rank MSE (low to high)
  mutate(Rank_MSE = rank(MSE, ties.method = "first")) %>%
  arrange(desc(DeltaK)) %>%  # Rank ΔK (high to low)
  mutate(Rank_DeltaK = rank(-DeltaK, ties.method = "first"))

# Reorder table by best average rank
ranking_table <- ranking_table %>%
  mutate(Average_Rank = (Rank_CE + Rank_MSE + Rank_DeltaK) / 3) %>%
  arrange(Average_Rank)

# Print ranked table
print(ranking_table) # 10, 9, 8, 7 (10 taxa)

##############################################################
#### Plotting 4 values of K in one figure ALternate one
##############################################################
##############################################################
# Load required libraries
library(ggplot2)
library(reshape2)
library(viridis)
library(gridExtra)
library(dplyr)
library(cowplot)  # for shared legend

# Define K values to analyze
K_values <- c(10, 9, 8, 7)

# Load metadata
metadata <- read.csv("affinis_metadata_floricomplex.csv", header = TRUE, sep = ",")

# Check that metadata contains required columns
if (!all(c("Sample_ID", "Taxon", "Population", "Order") %in% colnames(metadata))) {
  stop("Metadata must contain 'Sample_ID', 'Taxon', 'Population', and 'Order'.")
}

# Generate colorblind-friendly palette using viridis
generate_palette <- function(k) {
  viridis(n = k, option = "D")
}

# Function to generate sNMF barplot
plot_snmf <- function(k, show_legend = TRUE) {
  best_run <- which.min(cross.entropy(obj, K = k))
  qmatrix <- as.data.frame(Q(obj, K = k, run = best_run))
  
  # Rename columns for clarity
  colnames(qmatrix) <- paste0("Cluster_", seq_len(ncol(qmatrix)))
  
  # Add Sample_ID from metadata (assumes same row order)
  if (nrow(qmatrix) != nrow(metadata)) {
    stop("Mismatch between Q matrix rows and metadata rows.")
  }
  
  qmatrix$Sample_ID <- metadata$Sample_ID
  qmatrix$Population <- metadata$Population
  qmatrix$Taxon <- metadata$Taxon
  qmatrix$Order <- metadata$Order
  
  # Melt for ggplot
  q_long <- melt(qmatrix, id.vars = c("Sample_ID", "Taxon", "Population", "Order"),
                 variable.name = "Cluster", value.name = "Proportion")
  
  # Order samples
  q_long$Sample_ID <- factor(q_long$Sample_ID,
                             levels = metadata[order(metadata$Taxon, metadata$Population, metadata$Order), "Sample_ID"])
  
  # Colors
  k_colors <- generate_palette(k)
  
  # Build plot
  p <- ggplot(q_long, aes(x = Sample_ID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = k_colors) +
    labs(title = paste("sNMF Admixture Plot (K =", k, ")"),
         x = NULL, y = "Ancestry Proportion") +
    facet_grid(~ Taxon, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(1, "lines"),
      strip.text.x = element_text(size = 8),
      legend.position = if (show_legend) "bottom" else "none"
    )
  
  return(p)
}

# Create and save individual plots
plots_list <- list()
for (k in K_values) {
  p <- plot_snmf(k, show_legend = TRUE)
  plots_list[[as.character(k)]] <- p
  ggsave(filename = paste0("sNMF_K", k, "_plot.png"), plot = p, width = 14, height = 5, dpi = 300)
}

# Create versions without legends for composite
plots_nolegend <- lapply(K_values, function(k) plot_snmf(k, show_legend = FALSE))

# Extract one plot with legend
legend_plot <- plot_snmf(K_values[1], show_legend = TRUE)
shared_legend <- get_legend(legend_plot)

# Combine all plots into one figure with shared legend
composite_plot <- plot_grid(plotlist = plots_nolegend, ncol = 2, align = "v")
final_plot <- plot_grid(composite_plot, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot
# Save final composite figure
ggsave("sNMF_100reps_K7_K10_affinis_90.png", plot = final_plot, width = 16, height = 10, dpi = 300)


##############################################################
#### Plotting 4 values of K in one figure
##############################################################
##############################################################
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)
library(gridExtra)

# Define K values to analyze
K_values <- c(10, 9, 8, 7)  # Modify these K values as needed

# Load metadata
metadata <- read.csv("affinis_metadata_floricomplex.csv", header = TRUE, sep = ",")

# Ensure metadata contains necessary columns
if (!all(c("Sample_ID", "Taxon", "Population") %in% colnames(metadata))) {
  stop("Error: Metadata file must contain 'Sample_ID', 'Taxon', and 'Population' columns.")
}

# Extract sample names from the sNMF Q-matrix (first K value as reference)
qmatrix_samples <- rownames(Q(obj, K = K_values[1], run = 1))

# Identify missing samples (if any)
missing_from_Q <- setdiff(metadata$Sample_ID, qmatrix_samples)
missing_from_metadata <- setdiff(qmatrix_samples, metadata$Sample_ID)

# Print missing samples for debugging
if (length(missing_from_Q) > 0) {
  cat("⚠️ Samples in metadata but MISSING from Q-matrix:\n")
  print(missing_from_Q)
}

if (length(missing_from_metadata) > 0) {
  cat("⚠️ Samples in Q-matrix but MISSING from metadata:\n")
  print(missing_from_metadata)
}

# Function to generate a colorblind-friendly palette dynamically
generate_palette <- function(k) {
  viridis(n = k, option = "D")  # "D" is a colorblind-friendly viridis option
}

# Function to generate an sNMF bar plot for a given K
plot_snmf <- function(k) {
  best_run <- which.min(cross.entropy(obj, K = k))  # Select best run
  qmatrix <- as.data.frame(Q(obj, K = k, run = best_run))
  
  # Rename K cluster columns
  colnames(qmatrix) <- paste0("Cluster_", seq_len(ncol(qmatrix)))
  
  # Add metadata
  qmatrix$Sample_ID <- metadata$Sample_ID
  qmatrix$Population <- metadata$Population
  qmatrix$Taxon <- metadata$Taxon
  
  # Reshape for ggplot
  qmatrix_long <- melt(qmatrix, id.vars = c("Sample_ID", "Population", "Taxon"),
                       variable.name = "Cluster", value.name = "Proportion")
  
  # Ensure consistent sample order
  qmatrix_long$Sample_ID <- factor(qmatrix_long$Sample_ID, levels = metadata$Sample_ID)
  
  # Generate color palette
  k_colors <- generate_palette(k)
  
  # Create ggplot
  p <- ggplot(qmatrix_long, aes(x = Sample_ID, y = Proportion, fill = Cluster)) +  
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = k_colors) +  # Use dynamically generated colors
    theme_minimal() +
    labs(title = paste("sNMF Results for K =", k), fill = "Cluster") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~Taxon, scales = "free_x", nrow = 1)  # Group by Taxon
  
  return(p)
}

# Generate and save individual plots
for (k in K_values) {
  p <- plot_snmf(k)
  ###### Change name here for each different dataset for unique names. There is no overwrite of names : )
  ggsave(filename = paste0("sNMF_K_viridis_100reps_K220", k, ".png"), plot = p, width = 10, height = 4, dpi = 300)
  print(p)
}

# Generate combined figure with all K values
plots <- lapply(K_values, plot_snmf)
grid.arrange(grobs = plots, ncol = 2)  # Arrange in 2x2 grid

