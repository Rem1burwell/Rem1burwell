# =======================================================
# Genomic Comparison Pipeline: Fst, Gst, and Jost's D
# =======================================================
# Purpose: Compute and visualize pairwise population metrics
# using VCF data and metadata. Includes Mantel tests (IBD),
# regression summaries, and AMOVA analyses. Optimized for
# publication-ready figures and tables, with optional interactivity.

# ------------------------
# Load required libraries
# ------------------------
library(adegenet)
library(poppr)
library(pegas)
library(vegan)
library(dplyr)
library(future.apply)
library(ggplot2)
library(vcfR)
library(hierfstat)
library(doParallel)
library(reshape2)
library(gridExtra)
library(boot)
library(png)
library(grid)
library(mmod)        # Gst and Jost's D
library(viridis)     # Color palettes
library(geosphere)   # Geographic distance (IBD)
library(ggsignif)    # Significance stars
library(plotly)      # Interactive plots
library(sf)          # Simple Features for spatial geometry
library(ggspatial)   # For scale bar, north arrow

# ------------------------
# Parallel setup (optional)
# ------------------------
cores <- detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

setwd("~/Desktop/Affinis_genomics")
vcf <- read.vcfR("affinis_full_90.recode.vcf")  
# ------------------------
# Patch genind from vcfR2genind()
# ------------------------
patch_genind <- function(geno) {
  if (is.null(colnames(geno@tab))) {
    colnames(geno@tab) <- paste0("Locus_", seq_len(ncol(geno@tab)))
  }
  if (is.null(geno@loc.fac) || length(geno@loc.fac) != ncol(geno@tab)) {
    loci_names <- gsub("\\..*", "", colnames(geno@tab))
    geno@loc.fac <- factor(loci_names)
  }
  if (is.null(geno@loc.n.all) || length(geno@loc.n.all) != length(unique(geno@loc.fac))) {
    geno@loc.n.all <- as.integer(table(geno@loc.fac))
  }
  geno@all.names <- split(colnames(geno@tab), geno@loc.fac)
  if (length(unique(geno@ploidy)) > 1) {
    geno@ploidy <- rep(median(geno@ploidy), length(geno@ploidy))
  }
  return(geno)
}

# ------------------------
# Bootstrapped CI for regression
# ------------------------
boot_lm_ci <- function(model, data, indices) {
  d <- data[indices, ]
  coef(lm(formula(model), data = d))
}

# ------------------------
# Summary statistics table for regressions
# ------------------------
summarize_regression_stats <- function(results, output_file = "regression_summary.csv") {
  all_stats <- list()
  for (group in names(results)) {
    group_results <- results[[group]]$regressions
    if (is.null(group_results)) next
    for (model_name in names(group_results)) {
      model <- group_results[[model_name]]
      if (!inherits(model, "lm")) next
      stats <- summary(model)
      df <- data.frame(
        Grouping = group,
        Comparison = model_name,
        Intercept = coef(model)[1],
        Slope = coef(model)[2],
        R2 = stats$r.squared,
        Adj_R2 = stats$adj.r.squared,
        P_value = coef(stats)[2, 4]
      )
      all_stats[[paste(group, model_name, sep = "_")]] <- df
    }
  }
  final_df <- do.call(rbind, all_stats)
  write.csv(final_df, file = output_file, row.names = FALSE)
  return(final_df)
}

# ------------------------
# Render barplot of regression stats w/ significance stars
# ------------------------
render_regression_summary_plot <- function(summary_df, output_file = "regression_summary_plot.png") {
  plot_df <- summary_df %>%
    select(Grouping, Comparison, R2, Adj_R2, P_value) %>%
    pivot_longer(cols = c(R2, Adj_R2, P_value), names_to = "Metric", values_to = "Value")
  plot_df$Comparison <- gsub("_", " vs ", plot_df$Comparison)
  plot_df$stars <- cut(plot_df$P_value, 
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                       labels = c("***", "**", "*", ""))
  p <- ggplot(plot_df, aes(x = Comparison, y = Value, fill = Grouping)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    geom_text(aes(label = stars), vjust = -0.5, size = 4, position = position_dodge(width = 0.7)) +
    facet_wrap(~ Metric, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(title = "Regression Summary Statistics",
         x = "Comparison", y = "Value", fill = "Grouping") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(output_file, p, width = 12, height = 6)
  return(p)
}

# ------------------------
# Interactive regression summary (plotly)
# ------------------------
render_interactive_summary_plot <- function(summary_df) {
  library(plotly)
  plot_df <- summary_df %>%
    select(Grouping, Comparison, R2, Adj_R2, P_value) %>%
    pivot_longer(cols = c(R2, Adj_R2, P_value), names_to = "Metric", values_to = "Value")
  plot_df$Comparison <- gsub("_", " vs ", plot_df$Comparison)
  plot_ly(data = plot_df, x = ~Comparison, y = ~Value, color = ~Grouping, 
          type = 'bar', text = ~paste("Grouping:", Grouping, "<br>Value:", round(Value, 3)), 
          hoverinfo = 'text') %>%
    layout(barmode = 'group', 
           facet = ~Metric, 
           title = "Interactive Regression Summary Statistics",
           xaxis = list(title = "Comparison"),
           yaxis = list(title = "Value"))
}


# ------------------------
# Full analysis function: Compare Fst, Gst, and Jost's D with IBD and AMOVA
# ------------------------
# compareFstGstJostD_allComparisons <- function(vcf_file,
#                                               metadata,
#                                               group_cols = c("Taxon", "Population"),
#                                               output_dir = "Fst_Gst_results",
#                                               sig_threshold = 0.05,
#                                               boot_n = 100,
#                                               plot_sig_labels = TRUE,
#                                               do_amova = TRUE,
#                                               amova_dist_method = "nei",
#                                               min_group_size = 2,
#                                               max_missing_per_locus = 0.2,
#                                               max_missing_per_indiv = 0.2,
#                                               return_all = TRUE) {
#   if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#   
#   vcf <- read.vcfR(vcf_file)
#   geno <- vcfR2genind(vcf)
#   geno <- patch_genind(geno)
#   
#   # Filter missing data
#   missing_loci <- colMeans(is.na(geno@tab))
#   geno <- geno[, missing_loci <= max_missing_per_locus]
#   missing_indiv <- rowMeans(is.na(geno@tab))
#   geno <- geno[missing_indiv <= max_missing_per_indiv]
#   
#   metadata$SampleID <- as.character(metadata$Sample_ID)
#   gl_samples <- indNames(geno)
#   metadata <- metadata[metadata$SampleID %in% gl_samples, ]
#   metadata <- metadata[match(gl_samples, metadata$SampleID), ]
#   stopifnot(all(indNames(geno) %in% metadata$SampleID))
#   
#   # IBD support
#   has_coords <- all(c("Latitude", "Longitude") %in% colnames(metadata))
#   if (has_coords) {
#     coords <- metadata[, c("Longitude", "Latitude")]
#     geo_dist <- distm(coords, fun = distHaversine)
#     rownames(geo_dist) <- colnames(geo_dist) <- metadata$SampleID
#   }
#   
#   all_results <- list()
#   
#   for (group_col in group_cols) {
#     cat("
# Running comparison for:", group_col, "...
# ")
#     if (!group_col %in% colnames(metadata)) next
#     
#     pop_vector <- metadata[[group_col]]
#     strata_df <- data.frame(Pop = pop_vector)
#     rownames(strata_df) <- indNames(geno)
#     strata(geno) <- strata_df
#     pop(geno) <- factor(pop_vector)
#     geno <- patch_genind(geno)
#     
#     geno <- geno[pop(geno) %in% names(which(table(pop(geno)) >= min_group_size))]
#     if (length(unique(pop(geno))) < 2) next
#     
#     pairwise_fst <- hierfstat::pairwise.WCfst(geno)
#     pairwise_df <- melt(as.matrix(pairwise_fst), varnames = c("Pop1", "Pop2"), value.name = "Fst") %>%
#       filter(Pop1 != Pop2) %>%
#       mutate(across(c(Pop1, Pop2), as.character)) %>%
#       mutate(key = paste(pmin(Pop1, Pop2), pmax(Pop1, Pop2), sep = "_")) %>%
#       distinct(key, .keep_all = TRUE) %>%
#       select(-key)
#     
#     pairwise_gst <- tryCatch({
#       as.matrix(mmod::pairwise_Gst_Nei(geno))
#     }, error = function(e) matrix(NA, nrow = nrow(pairwise_fst), ncol = ncol(pairwise_fst)))
#     
#     pairwise_jostd <- tryCatch({
#       as.matrix(mmod::pairwise_D(geno))
#     }, error = function(e) matrix(NA, nrow = nrow(pairwise_fst), ncol = ncol(pairwise_fst)))
#     
#     colnames(pairwise_gst) <- rownames(pairwise_gst) <- colnames(pairwise_jostd) <- rownames(pairwise_jostd) <- rownames(pairwise_fst)
#     
#     melt_gst <- melt(pairwise_gst, varnames = c("Pop1", "Pop2"), value.name = "Gst")
#     melt_jost <- melt(pairwise_jostd, varnames = c("Pop1", "Pop2"), value.name = "JostD")
#     
#     merged_all <- Reduce(function(x, y) merge(x, y, by = c("Pop1", "Pop2")), list(pairwise_df, melt_gst, melt_jost))
#     
#     if (has_coords) {
#       merged_all$GeoDist <- mapply(function(x, y) {
#         i <- which(metadata$SampleID == x)
#         j <- which(metadata$SampleID == y)
#         if (length(i) > 0 && length(j) > 0) geo_dist[i, j] else NA
#       }, merged_all$Pop1, merged_all$Pop2)
#     }
#     
#     reg_models <- list()
#     reg_ci <- list()
#     regressions <- list(
#       Fst_Gst = list(x = "Fst", y = "Gst"),
#       Fst_Jost = list(x = "Fst", y = "JostD"),
#       Gst_Jost = list(x = "Gst", y = "JostD")
#     )
#     if (has_coords) {
#       regressions$IBD_Fst <- list(x = "GeoDist", y = "Fst")
#       regressions$IBD_JostD <- list(x = "GeoDist", y = "JostD")
#     }
#     
#     for (name in names(regressions)) {
#       x <- regressions[[name]]$x
#       y <- regressions[[name]]$y
#       subdata <- merged_all %>% filter(!is.na(.data[[x]]) & !is.na(.data[[y]]))
#       if (nrow(subdata) < 2) next
#       model <- lm(as.formula(paste(y, "~", x)), data = subdata)
#       reg_models[[name]] <- model
#       ci <- tryCatch(boot::boot(data = subdata, statistic = boot_lm_ci, R = boot_n, model = model), error = function(e) NULL)
#       reg_ci[[name]] <- ci
#       plot <- ggplot(subdata, aes_string(x = x, y = y, color = "Pop1", shape = "Pop1")) +
#         geom_point(size = 3) +
#         geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
#         scale_color_viridis_d(option = "D", end = 0.85) +
#         scale_shape_manual(values = 1:length(unique(subdata$Pop1))) +
#         labs(title = paste(group_col, ":", x, "vs", y), x = x, y = y, color = "Group", shape = "Group") +
#         theme_minimal()
#       ggsave(file.path(output_dir, paste0("reg_", tolower(name), "_", group_col, ".png")), plot, width = 8, height = 6)
#     }
#     
#     amova_res <- NULL
#     if (do_amova) {
#       distmat <- tryCatch({
#         if (amova_dist_method == "nei") poppr::bruvo.dist(geno)
#         else dist(tab(geno, NA.method = "mean"))
#       }, error = function(e) NULL)
#       if (!is.null(distmat)) {
#         amova_df <- data.frame(Pop = pop(geno))
#         if (!is.factor(amova_df$Pop)) amova_df$Pop <- as.factor(amova_df$Pop)
#         try({
#           amova_res <- pegas::amova(distmat ~ Pop, data = amova_df, nperm = boot_n)
#           sink(file.path(output_dir, paste0("amova_", group_col, ".txt")))
#           print(summary(amova_res))
#           sink()
#         }, silent = TRUE)
#       }
#     }
#     
#     write.csv(merged_all, file = file.path(output_dir, paste0("pairwise_metrics_", group_col, ".csv")), row.names = FALSE)
#     all_results[[group_col]] <- list(
#       summary = merged_all,
#       regressions = reg_models,
#       amova = amova_res,
#       ci = reg_ci
#     )
#   }
#   if (return_all) return(all_results) else invisible(NULL)
# }
patch_genind <- function(geno) {
  if (is.null(colnames(geno@tab))) {
    colnames(geno@tab) <- paste0("Locus_", seq_len(ncol(geno@tab)))
  }
  if (is.null(geno@loc.fac) || length(geno@loc.fac) != ncol(geno@tab)) {
    loci_names <- gsub("\\..*", "", colnames(geno@tab))
    geno@loc.fac <- factor(loci_names)
  }
  if (is.null(geno@loc.n.all) || length(geno@loc.n.all) != length(unique(geno@loc.fac))) {
    geno@loc.n.all <- as.integer(table(geno@loc.fac))
  }
  geno@all.names <- split(colnames(geno@tab), geno@loc.fac)
  if (length(unique(geno@ploidy)) > 1) {
    geno@ploidy <- rep(median(geno@ploidy), length(geno@ploidy))
  }
  return(geno)
}

boot_lm_ci <- function(model, data, indices) {
  d <- data[indices, ]
  coef(lm(formula(model), data = d))
}

summarize_regression_stats <- function(results, output_file = "regression_summary.csv") {
  all_stats <- list()
  
  for (group in names(results)) {
    group_results <- results[[group]]$regressions
    if (is.null(group_results)) next
    
    for (model_name in names(group_results)) {
      model <- group_results[[model_name]]
      if (!inherits(model, "lm")) next
      
      stats <- summary(model)
      df <- data.frame(
        Grouping = group,
        Comparison = model_name,
        Intercept = coef(model)[1],
        Slope = coef(model)[2],
        R2 = stats$r.squared,
        Adj_R2 = stats$adj.r.squared,
        P_value = coef(stats)[2, 4]
      )
      all_stats[[paste(group, model_name, sep = "_")]] <- df
    }
  }
  
  final_df <- do.call(rbind, all_stats)
  write.csv(final_df, file = output_file, row.names = FALSE)
  return(final_df)
}

compareFstGstJostD_allComparisons <- function(vcf_file, metadata, group_cols = c("Taxon", "Population"), output_dir = "Fst_Gst_results", sig_threshold = 0.05, boot_n = 100, plot_sig_labels = TRUE, do_amova = TRUE, amova_dist_method = "nei", min_group_size = 2, max_missing_per_locus = 0.2, max_missing_per_indiv = 0.2, return_all = TRUE) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  vcf <- read.vcfR(vcf_file)
  geno <- vcfR2genind(vcf)
  geno <- patch_genind(geno)
  
  missing_loci <- colMeans(is.na(geno@tab))
  geno <- geno[, missing_loci <= max_missing_per_locus]
  missing_indiv <- rowMeans(is.na(geno@tab))
  geno <- geno[missing_indiv <= max_missing_per_indiv]
  
  metadata$SampleID <- as.character(metadata$Sample_ID)
  gl_samples <- indNames(geno)
  metadata <- metadata[metadata$SampleID %in% gl_samples, ]
  metadata <- metadata[match(gl_samples, metadata$SampleID), ]
  stopifnot(all(indNames(geno) %in% metadata$SampleID))
  
  has_coords <- all(c("Latitude", "Longitude") %in% colnames(metadata))
  if (has_coords) {
    coord_metadata <- metadata[complete.cases(metadata[, c("Latitude", "Longitude")]), ]
    coords <- coord_metadata[, c("Longitude", "Latitude")]
    geo_dist <- distm(coords, fun = distHaversine)
    rownames(geo_dist) <- colnames(geo_dist) <- coord_metadata$SampleID
  }
  
  all_results <- list()
  
  for (group_col in group_cols) {
    cat("\nRunning comparison for:", group_col, "...\n")
    if (!group_col %in% colnames(metadata)) next
    
    pop_vector <- metadata[[group_col]]
    strata_df <- data.frame(Pop = pop_vector)
    rownames(strata_df) <- indNames(geno)
    strata(geno) <- strata_df
    pop(geno) <- factor(pop_vector)
    geno <- patch_genind(geno)
    
    geno <- geno[pop(geno) %in% names(which(table(pop(geno)) >= min_group_size))]
    if (length(unique(pop(geno))) < 2) next
    
    pairwise_fst <- hierfstat::pairwise.WCfst(geno)
    pairwise_df <- melt(as.matrix(pairwise_fst), varnames = c("Pop1", "Pop2"), value.name = "Fst") %>%
      filter(Pop1 != Pop2) %>%
      mutate(across(c(Pop1, Pop2), as.character)) %>%
      mutate(key = paste(pmin(Pop1, Pop2), pmax(Pop1, Pop2), sep = "_")) %>%
      distinct(key, .keep_all = TRUE) %>%
      select(-key)
    
    pairwise_gst <- tryCatch({
      as.matrix(mmod::pairwise_Gst_Nei(geno))
    }, error = function(e) matrix(NA, nrow = nrow(pairwise_fst), ncol = ncol(pairwise_fst)))
    
    pairwise_jostd <- tryCatch({
      as.matrix(mmod::pairwise_D(geno))
    }, error = function(e) matrix(NA, nrow = nrow(pairwise_fst), ncol = ncol(pairwise_fst)))
    
    colnames(pairwise_gst) <- rownames(pairwise_gst) <- colnames(pairwise_jostd) <- rownames(pairwise_jostd) <- rownames(pairwise_fst)
    
    melt_gst <- melt(pairwise_gst, varnames = c("Pop1", "Pop2"), value.name = "Gst")
    melt_jost <- melt(pairwise_jostd, varnames = c("Pop1", "Pop2"), value.name = "JostD")
    
    merged_all <- Reduce(function(x, y) merge(x, y, by = c("Pop1", "Pop2")), list(pairwise_df, melt_gst, melt_jost))
    
    if (has_coords) {
      merged_all$GeoDist <- mapply(function(x, y) {
        if (x %in% rownames(geo_dist) && y %in% colnames(geo_dist)) {
          geo_dist[x, y]
        } else {
          NA
        }
      }, merged_all$Pop1, merged_all$Pop2)
    }
    
    reg_models <- list()
    reg_ci <- list()
    regressions <- list(
      Fst_Gst = list(x = "Fst", y = "Gst"),
      Fst_Jost = list(x = "Fst", y = "JostD"),
      Gst_Jost = list(x = "Gst", y = "JostD")
    )
    if (has_coords) {
      regressions$IBD_Fst <- list(x = "GeoDist", y = "Fst")
      regressions$IBD_JostD <- list(x = "GeoDist", y = "JostD")
    }
    
    for (name in names(regressions)) {
      x <- regressions[[name]]$x
      y <- regressions[[name]]$y
      subdata <- merged_all %>% filter(!is.na(.data[[x]]) & !is.na(.data[[y]]))
      if (nrow(subdata) < 2) next
      model <- lm(as.formula(paste(y, "~", x)), data = subdata)
      reg_models[[name]] <- model
      ci <- tryCatch(boot::boot(data = subdata, statistic = boot_lm_ci, R = boot_n, model = model), error = function(e) NULL)
      reg_ci[[name]] <- ci
      plot <- ggplot(subdata, aes_string(x = x, y = y, color = "Pop1", shape = "Pop1")) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
        scale_color_viridis_d(option = "D", end = 0.85) +
        scale_shape_manual(values = 1:length(unique(subdata$Pop1))) +
        labs(title = paste(group_col, ":", x, "vs", y), x = x, y = y, color = "Group", shape = "Group") +
        theme_minimal()
      ggsave(file.path(output_dir, paste0("reg_", tolower(name), "_", group_col, ".png")), plot, width = 8, height = 6)
    }
    
    amova_res <- NULL
    if (do_amova) {
      distmat <- tryCatch({
        if (amova_dist_method == "nei") poppr::bruvo.dist(geno)
        else dist(tab(geno, NA.method = "mean"))
      }, error = function(e) NULL)
      if (!is.null(distmat)) {
        amova_df <- data.frame(Pop = pop(geno))
        if (!is.factor(amova_df$Pop)) amova_df$Pop <- as.factor(amova_df$Pop)
        try({
          amova_res <- pegas::amova(distmat ~ Pop, data = amova_df, nperm = boot_n)
          sink(file.path(output_dir, paste0("amova_", group_col, ".txt")))
          print(summary(amova_res))
          sink()
        }, silent = TRUE)
      }
    }
    
    write.csv(merged_all, file = file.path(output_dir, paste0("pairwise_metrics_", group_col, ".csv")), row.names = FALSE)
    all_results[[group_col]] <- list(
      summary = merged_all,
      regressions = reg_models,
      amova = amova_res,
      ci = reg_ci
    )
  }
  
  if (return_all) return(all_results) else invisible(NULL)
}




# ------------------------
# Run full pipeline
# ------------------------
results <- compareFstGstJostD_allComparisons(
  vcf_file = "affinis_full_90.recode.vcf",
  metadata = read.csv("affinis_metadata_floricomplex.csv"),
  group_cols = c("Taxon", "Population"),
  output_dir = "Fst_Gst_comparisons_Affinis90",
  sig_threshold = 0.05,
  boot_n = 100,
  plot_sig_labels = TRUE,
  do_amova = TRUE,
  amova_dist_method = "nei",
  max_missing_per_locus = 0.05,
  max_missing_per_indiv = 0.10,
  min_group_size = 2,
  return_all = TRUE
)

# =======================================================
# Genomic Comparison Pipeline: Fst, Gst, and Jost's D
# =======================================================
# Purpose: Compute and visualize pairwise population metrics
# using VCF data and metadata. Includes Mantel tests (IBD),
# regression summaries, and AMOVA analyses. Optimized for
# publication-ready figures and tables, with optional interactivity.

# ------------------------
# Function to render plots from results (selectable)
# ------------------------
render_all_summary_plots <- function(results, metadata, output_dir = "Fst_Gst_comparisons_Affinis90", which_plots = c("regressions", "mantel", "map")) {
  if ("regressions" %in% which_plots) {
    summary_df <- summarize_regression_stats(results)
    plot_paths <- list.files(output_dir, pattern = "^reg_.*\\.png$", full.names = TRUE)
    plot_grobs <- lapply(plot_paths, function(path) {
      tryCatch({ rasterGrob(png::readPNG(path), interpolate = TRUE) }, error = function(e) NULL)
    })
    plot_grobs <- Filter(Negate(is.null), plot_grobs)
    if (length(plot_grobs) > 0) gridExtra::grid.arrange(grobs = plot_grobs, ncol = 2)
  }
  if ("mantel" %in% which_plots) {
    mantel_plot <- render_mantel_table_plot(results)
    print(mantel_plot)
  }
  if ("map" %in% which_plots) {
    map <- map_pairwise_metrics(results, metadata, metric = "Fst", output_file = file.path(output_dir, "pairwise_map_plot.png"))
    print(map)
  }
}


# Render only regression plots
render_all_summary_plots(results, metadata, which_plots = c("regressions"))

# Render only the map
render_all_summary_plots(results, metadata, which_plots = c("map"))

# Render all plots
render_all_summary_plots(results, metadata, which_plots = c("regressions", "mantel", "map"))


# ------------------------
# View and save summaries
# ------------------------
summary_df <- summarize_regression_stats(results)
print(summary_df)
render_regression_summary_plot(summary_df)
render_interactive_summary_plot(summary_df)

# ------------------------
# Display plots as grobs
# ------------------------
plot_paths <- list.files("Fst_Gst_comparisons_Affinis90", pattern = "^reg_.*\\.png$", full.names = TRUE)
plot_grobs <- lapply(plot_paths, function(path) rasterGrob(readPNG(path), interpolate = TRUE))
grid.arrange(grobs = plot_grobs, ncol = 2)

# ------------------------
# Ploty interative plotting
# ------------------------
render_interactive_summary_plot(summary_df)

# ------------------------
# Mapping function: pairwise Fst values with spatial lines
# ------------------------
map_pairwise_metrics <- function(results, metadata, metric = "Fst", output_file = "pairwise_map_plot.png") {
  require(ggplot2)
  require(sf)
  require(ggspatial)
  
  coords <- metadata[, c("SampleID", "Longitude", "Latitude")]
  coords <- coords[complete.cases(coords), ]
  coords_sf <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = 4326)
  
  lines <- list()
  
  for (group in names(results)) {
    df <- results[[group]]$summary
    if (!(metric %in% colnames(df))) next
    df <- df[complete.cases(df[, c("Pop1", "Pop2", metric)]), ]
    
    for (i in seq_len(nrow(df))) {
      p1 <- coords_sf[coords_sf$SampleID == df$Pop1[i], ]
      p2 <- coords_sf[coords_sf$SampleID == df$Pop2[i], ]
      if (nrow(p1) == 1 && nrow(p2) == 1) {
        line <- st_sfc(st_linestring(rbind(st_coordinates(p1), st_coordinates(p2))), crs = 4326)
        lines[[length(lines) + 1]] <- st_sf(
          geometry = line,
          Grouping = group,
          Value = df[[metric]][i]
        )
      }
    }
  }
  
  if (length(lines) == 0) return(NULL)
  
  all_lines <- do.call(rbind, lines)
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  p <- ggplot() +
    geom_sf(data = world, fill = "gray95", color = "gray80") +
    geom_sf(data = all_lines, aes(color = Value), size = 1) +
    scale_color_viridis_c(option = "D") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Pairwise", metric, "Relationships"), color = metric) +
    annotation_scale(location = "bl") +
    annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_minimal())
  
  ggsave(output_file, p, width = 10, height = 6)
  return(p)
}

map <- map_pairwise_metrics(results, metadata, metric = "Fst", output_file = "pairwise_map_plot.png")

# ------------------------
# Modular summary plot rendering
# ------------------------
render_all_summary_plots <- function(results, metadata, output_dir = "Fst_Gst_comparisons_Affinis90", which_plots = c("regressions", "mantel", "map", "amova"), metric = "Fst", interactive = FALSE) {
  if ("regressions" %in% which_plots) {
    summary_df <- summarize_regression_stats(results)
    if (interactive) {
      print(render_interactive_summary_plot(summary_df))
    } else {
      print(render_regression_summary_plot(summary_df))
    }
    plot_paths <- list.files(output_dir, pattern = "^reg_.*\\.png$", full.names = TRUE)
    plot_grobs <- lapply(plot_paths, function(path) {
      tryCatch({ rasterGrob(png::readPNG(path), interpolate = TRUE) }, error = function(e) NULL)
    })
    plot_grobs <- Filter(Negate(is.null), plot_grobs)
    if (length(plot_grobs) > 0) gridExtra::grid.arrange(grobs = plot_grobs, ncol = 2)
  }
  
  if ("mantel" %in% which_plots) {
    mantel_plot <- render_mantel_table_plot(results)
    if (!is.null(mantel_plot)) print(mantel_plot)
  }
  
  if ("map" %in% which_plots) {
    if (!all(c("SampleID", "Longitude", "Latitude") %in% colnames(metadata))) {
      warning("Cannot render map: metadata is missing required columns (SampleID, Longitude, Latitude).")
    } else {
      map <- tryCatch({
        map_pairwise_metrics(results, metadata, metric = metric, output_file = file.path(output_dir, paste0("pairwise_map_plot_", metric, ".png")))
      }, error = function(e) {
        warning("Map rendering failed: ", conditionMessage(e))
        return(NULL)
      })
      if (!is.null(map)) print(map)
    }
  }
  
  if ("amova" %in% which_plots) {
    amova_files <- list.files(output_dir, pattern = "^amova_.*\\.txt$", full.names = TRUE)
    if (length(amova_files) == 0) {
      message("No AMOVA summaries found.")
    } else {
      for (f in amova_files) {
        cat("\n--- AMOVA Summary:", basename(f), "---\n")
        cat(readLines(f), sep = "\n")
      }
    }
  }
}

##### Example usgae

render_all_summary_plots(results, metadata, which_plots = c("regressions", "mantel", "map", "amova"))

render_all_summary_plots(results, metadata, which_plots = c("map"), metric = "JostD")

# End of script
