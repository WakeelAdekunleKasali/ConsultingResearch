daat <- readRDS("/Users/wakeelkasali/Desktop/PROJECT_STAT_540/Overall/final_seurat_list.rds")

names(daat)

ctrl0 <- daat[["D0_Control"]]
asd0  <- daat[["D0_ASD"]]
ctrl6 <- daat[["D6_Control"]]
asd6  <- daat[["D6_ASD"]]
ctrl30 <- daat[["D30_Control"]]
asd30  <- daat[["D30_ASD"]]
ctrl60 <- daat[["D60_Control"]]
asd60  <- daat[["D60_ASD"]]
ctrl100 <- daat[["D100_Control"]]
asd100  <- daat[["D100_ASD"]]


# Function to compute E/I ratio from a Seurat object
compute_ei_ratio <- function(seurat_obj,
                             exc_genes,
                             inh_genes,
                             assay = "SCT", layer = "data") {
  # Set default assay
  DefaultAssay(seurat_obj) <- assay
  
  # Find which genes are present in this object
  present_exc <- intersect(exc_genes,
                           rownames(seurat_obj))
  present_inh <- intersect(inh_genes,
                           rownames(seurat_obj))
  
  # If none of the markers are present, return NA
  if (length(present_exc) == 0 & length(present_inh) == 0) {
    warning("No excitatory or inhibitory genes found in object.")
    return(NA)
  }
  
  # Fetch expression matrix from the layer
  expr_data <- FetchData(seurat_obj,
                         vars = unique(c(present_exc, present_inh)),
                         layer = layer)
  
  # Classify cell types
  seurat_obj$Etype <- ifelse(rowSums(expr_data[, present_exc, drop = FALSE]) > 0,
                             "Excitatory",
                             ifelse(rowSums(expr_data[, present_inh, drop = FALSE]) > 0,
                                    "Inhibitory", "Unclassified"))
  
  # Count cell types
  etype_counts <- table(seurat_obj$Etype)
  
  # Extract counts
  n_exc <- ifelse("Excitatory" %in% names(etype_counts),
                  etype_counts["Excitatory"], NA)
  n_inh <- ifelse("Inhibitory" %in% names(etype_counts),
                  etype_counts["Inhibitory"], NA)
  
  # Compute E/I ratio
  ratio <- as.numeric(n_exc) / as.numeric(n_inh)
  return(ratio)


  # Define gene lists globally (so we don't repeat them for every call)
  exc_genes <- c("SLC17A6", "SLC17A7", "TBR1", "BCL11B", "SATB2", "FEZF2")
  inh_genes <- c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "LHX6")
  
  # Create a dataframe with E/I ratios
  ei_df <- data.frame(
    Timepoint = c("D30", "D30", "D60", "D60", "D100", "D100"),
    Group = c("Control", "ASD", "Control", "ASD", "Control", "ASD"),
    EI_Ratio = c(
      compute_ei_ratio(ctrl30, exc_genes, inh_genes),
      compute_ei_ratio(asd30, exc_genes, inh_genes),
      compute_ei_ratio(ctrl60, exc_genes, inh_genes),
      compute_ei_ratio(asd60, exc_genes, inh_genes),
      compute_ei_ratio(ctrl100, exc_genes, inh_genes),
      compute_ei_ratio(asd100, exc_genes, inh_genes)
    )
  )
  
  ei_df$Timepoint <- factor(ei_df$Timepoint, levels = c("D30", "D60", "D100"))
  
  ggplot(ei_df, aes(x = Timepoint, y = EI_Ratio, color = Group, group = Group)) +
    geom_line(size = 1.0) +
    geom_point(size = 1.5) +
    labs(
      title = "E/I Ratio Across Developmental Timepoints",
      x = "Developmental Timepoint",
      y = "Excitatory / Inhibitory Ratio"
    ) +
    theme_bw() +
    scale_color_manual(values = c("Control" = "#1f77b4", "ASD" = "#d62728")) +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  # Define gene lists (if not defined earlier)
  exc_genes <- c("SLC17A6", "SLC17A7", "TBR1", "BCL11B", "SATB2", "FEZF2")
  inh_genes <- c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "LHX6")
  
  # Step 1: Build dataframe dynamically using the compute_ei_ratio function
  ei_data <- tibble::tibble(
    Timepoint = rep(c(30, 60, 100), each = 2),
    Group = rep(c("Control", "ASD"), times = 3),
    EIRatio = c(
      compute_ei_ratio(ctrl30, exc_genes, inh_genes),
      compute_ei_ratio(asd30, exc_genes, inh_genes),
      compute_ei_ratio(ctrl60, exc_genes, inh_genes),
      compute_ei_ratio(asd60, exc_genes, inh_genes),
      compute_ei_ratio(ctrl100, exc_genes, inh_genes),
      compute_ei_ratio(asd100, exc_genes, inh_genes)
    )
  )
  
  
  ei_data_wide <- ei_data %>%
    pivot_wider(names_from = Group, values_from = EIRatio) %>%
    mutate(Diff = ASD - Control)
  
  # Step 3: Plot divergence
  ggplot(ei_data_wide, aes(x = factor(Timepoint, levels = c(30, 60, 100)), y = Diff)) +
    geom_point(size = 2.5, color = "blue") +
    geom_hline(yintercept = 0, linewidth = 1.2, linewidth = 2) +
    ylab("E/I Ratio Difference (ASD - Control)") +
    xlab("Developmental Timepoint") +
    ggtitle("Divergence in E/I Ratio Across Time") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )