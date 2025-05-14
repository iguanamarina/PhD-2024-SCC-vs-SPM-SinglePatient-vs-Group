############################## ################### ############## ### 
##
## Script name: 1. MASTERSCRIPT FOR SINGLE PATIENTS VS GROUP.r
##
## Purpose of script: Estimate sensitivity, specificity, and predictive
##                    metrics for SCCs and SPM in a 1 vs Group setup using
##                    simulated PET images and known ROIs.
##
##                    In this setting, a single simulated Alzheimer patient
##                    is compared against a control group. To enable SCC-based
##                    inference, artificial variability is added to the single
##                    subject using Poisson clones via generatePoissonClones().
##
##                    True regions are known, allowing performance comparison
##                    with SPM binary detection files. This script follows the
##                    same logic and structure as the group vs group script.
##
## Notes: Some computations (e.g., SCC estimation) can be skipped by using
##        precomputed `.RData` files available in the z35/ subfolder.
##
## Date Created: 2022-01-10
## Date Updated: 2025-05-12
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://juan-arias.xyz
##   
############################## ################### ############## ### 

### ==================================================== ###
###                     1) PREAMBLE                      ###
### ==================================================== ###

#* Set working directory
setwd("~/GitHub/PhD-2024-SCC-vs-SPM-SinglePatient-vs-Group")

#* Options
options(scipen = 6, digits = 4)

#* Define axial slice to analyze
paramZ <- 35

#* Define CRAN packages
cran_packages <- c(
  "ggplot2", "patchwork", "fields", "viridis",     # Plotting
  "knitr", "kableExtra", "tibble", "magrittr",     # Tables & reports
  "tidyr", "dplyr", "scales",                      # Data wrangling
  "xtable", "neuroSCC", "threadr", "stringr",      # Core and utility packages
  "readr"                                          # File handling
)

#* Load CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#* Load drat-hosted packages (non-CRAN)
drat_packages <- c("Triangulation", "ImageSCC", "BPST")
for (pkg in drat_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = c(neuroSCC = "https://iguanamarina.github.io/neuroSCC-drat",
                CRAN = "https://cloud.r-project.org"),
      type = "source"
    )
  }
  library(pkg, character.only = TRUE)
}

#* Clean up
rm(pkg, cran_packages, drat_packages)

### ==================================================== ###
###             2) CONTOURS OF NEURO-DATA               ###
### ==================================================== ###

#* Define file paths
mask_path <- "Auxiliary Files/new_mask.nii"
triangulation_path <- file.path("Results", "z35", "contour35.RData")

#* Only run if triangulation file does not exist
if (!file.exists(triangulation_path)) {
  
  message("[INFO] Generating triangulation from mask...")
  
  # Extract outer brain contour at level 0
  contours <- neuroSCC::neuroContour(mask_path, levels = 0)
  triangulation <- Triangulation::TriMesh(contours[[1]], n = 15)
  
  save(triangulation, file = triangulation_path)
  message("[INFO] Triangulation saved to: ", triangulation_path)
  
} else {
  message("[INFO] Triangulation already exists. Skipping generation.")
}

#* Clean up
rm(mask_path, triangulation_path)

### ==================================================== ###
###     3) CREATE SCC MATRIX FOR CONTROL GROUP          ###
### ==================================================== ###

#* NOTE: This matrix is assumed to have been generated in the simulations script
#        (MASTERSCRIPT for group comparisons) and saved to Results/z35/SCC_CN.RData.
#        If not already available, uncomment the code below to recompute it.

#* Load precomputed matrix
matrixCN_path <- file.path("Results", "z35", "SCC_CN.RData")

if (file.exists(matrixCN_path)) {
  message("[INFO] Loading SCC matrix for Control group...")
  load(matrixCN_path)  # loads: SCC_CN
} else {
  stop("[ERROR] SCC_CN.RData not found in Results/z35/. Please generate it using the group vs group script.")
}

# #* OPTIONAL: Recompute matrix from scratch if needed
# # -----------------------------------------------------
# message("[INFO] Creating SCC matrix for Control group from scratch...")
# 
# setwd("PETimg_masked for simulations")
# 
# databaseCN <- neuroSCC::databaseCreator(
#   pattern = "_w00_",
#   control = TRUE,
#   quiet = FALSE
# )
# 
# setwd("..")
# 
# SCC_CN <- neuroSCC::matrixCreator(
#   database = databaseCN,
#   paramZ = paramZ,
#   quiet = FALSE
# )
# 
# SCC_CN <- neuroSCC::meanNormalization(SCC_CN, quiet = FALSE)
# 
# save(SCC_CN, file = matrixCN_path)
# message("[INFO] Control matrix saved to: ", matrixCN_path)

#* Clean up
rm(matrixCN_path)


### ==================================================== ###
###       4) TRIANGULATION AND SCC PARAMETERS           ###
### ==================================================== ###

#* Load triangulation from Results/z35
load(file.path("Results", "z35", "contour35.RData"))  # loads: triangulation

#* Extract triangulation components
Brain.V <- VT[[1]]
Brain.Tr <- VT[[2]]

V.est  <- as.matrix(Brain.V)
Tr.est <- as.matrix(Brain.Tr)
V.band <- V.est
Tr.band <- Tr.est

#* Define regions and ROIs for iteration
regions <- c("roiAD", "w32", "w79", "w214", "w271", "w413")
rois    <- c(1, 2, 4, 6, 8)

#* Define SCC hyperparameters (Wang et al.)
d.est      <- 5
d.band     <- 2
r          <- 1
lambda     <- 10^{seq(-6, 3, 0.5)}
alpha.grid <- c(0.10, 0.05, 0.01)


### ==================================================== ###
###     5) SCC EVALUATION (1 vs Group)                  ###
### ==================================================== ###

#* Paths
spmDir     <- file.path("Results", "z35", "1vsGroup", "SPM")
resultsDir <- file.path("Results", "z35", "results")
roiDir     <- "roisNormalizadas"
maskPath   <- "Auxiliary Files/new_mask.nii"
dims       <- neuroSCC::getDimensions(maskPath)

#* List SPM files
spmFiles <- list.files(spmDir, pattern = "^binary_swwwC\\d+_\\w+_[148]\\.nii$", full.names = TRUE)

#* Extract metadata
meta <- stringr::str_match(basename(spmFiles),
                           pattern = "^binary_swww(C\\d+)_([a-zA-Z0-9]+)_([148])\\.nii$")
colnames(meta) <- c("full", "subjectID", "region", "roi")
meta <- as.data.frame(meta[, -1], stringsAsFactors = FALSE)
meta$roi <- as.integer(meta$roi)

#* Filter target regions and ROIs
meta <- meta |>
  dplyr::filter(region %in% c("w32", "w214", "w271", "roiAD"),
                roi %in% c(1, 4, 8))

#* Loop through region × roi
for (regionName in unique(meta$region)) {
  for (roiLevel in unique(meta$roi)) {
    
    message("[INFO] Evaluating SPM for region ", regionName, ", ROI ", roiLevel, "...")
    
    # Subset all matching subject files
    subsetMeta <- meta |> dplyr::filter(region == regionName, roi == roiLevel)
    if (nrow(subsetMeta) == 0) next
    
    outputDir  <- file.path(resultsDir, paste0("ROI", roiLevel))
    outputCSV  <- file.path(outputDir, paste0("sens_esp_SPM_", regionName, "_", roiLevel, ".csv"))
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    
    metricsList <- list()
    
    for (k in seq_len(nrow(subsetMeta))) {
      row <- subsetMeta[k, ]
      
      niftiFile     <- file.path(spmDir, paste0("binary_swww", row$subjectID, "_", row$region, "_", row$roi, ".nii"))
      roiMaskRegion <- ifelse(row$region == "roiAD", "wroiAD", row$region)
      roiFile       <- file.path(roiDir, paste0("wwwx", roiMaskRegion,
                                                "_redim_crop_squ_flipLR_newDim_",
                                                row$subjectID, ".nii"))
      
      if (!file.exists(niftiFile)) {
        warning("[WARNING] SPM file missing: ", niftiFile)
        next
      }
      if (!file.exists(roiFile)) {
        warning("[WARNING] ROI file missing: ", roiFile)
        next
      }
      
      # Load SPM and ROI points
      detectedPoints <- neuroSCC::getSPMbinary(niftiFile, paramZ = paramZ)
      truePoints <- neuroSCC::processROIs(
        roiFile = roiFile,
        region  = roiMaskRegion,
        number  = row$subjectID,
        save    = FALSE,
        verbose = FALSE
      )
      
      # Skip if no signal at this slice
      sliceTruePoints <- subset(truePoints, z == paramZ & pet == 1)
      if (nrow(sliceTruePoints) == 0) {
        message("[WARNING] No ROI voxels at z = ", paramZ, " for ", row$subjectID, " — skipping.")
        next
      }
      
      # Evaluate metrics
      subjectMetrics <- neuroSCC::calculateMetrics(
        detectedPoints = detectedPoints,
        truePoints     = truePoints,
        totalCoords    = dims,
        regionName     = paste0(row$region, "_", row$roi)
      )
      
      subjectMetrics$subject <- row$subjectID
      metricsList[[length(metricsList) + 1]] <- subjectMetrics
      
      rm(detectedPoints, truePoints, sliceTruePoints, subjectMetrics)
    }
    
    # Write full batch of results (guaranteed consistent structure)
    if (length(metricsList) > 0) {
      allMetrics <- do.call(rbind, metricsList)
      readr::write_csv(allMetrics, outputCSV)
      message("[INFO] Metrics written to: ", outputCSV)
    } else {
      warning("[WARNING] No SPM metrics computed for ", regionName, " ROI ", roiLevel)
    }
    
    rm(metricsList, allMetrics)
  }
}

#* Clean up
rm(spmDir, resultsDir, roiDir, maskPath, dims, spmFiles, meta,
   outputDir, outputCSV, roiFile, niftiFile, roiMaskRegion)


### ==================================================== ###
###     6) SPM EVALUATION (1 vs Group)                 ###
### ==================================================== ###

#* Paths
spmDir     <- file.path("Results", "z35", "1vsGroup", "SPM")
resultsDir <- file.path("Results", "z35", "results")
roiDir     <- "roisNormalizadas"
maskPath   <- "Auxiliary Files/new_mask.nii"
dims       <- neuroSCC::getDimensions(maskPath)

#* List available SPM files
spmFiles <- list.files(spmDir, pattern = "^binary_swwwC\\d+_\\w+_[148]\\.nii$", full.names = TRUE)

#* Extract metadata
meta <- stringr::str_match(basename(spmFiles),
                           pattern = "^binary_swww(C\\d+)_([a-zA-Z0-9]+)_([148])\\.nii$")
colnames(meta) <- c("full", "subjectID", "region", "roi")
meta <- as.data.frame(meta[, -1], stringsAsFactors = FALSE)
meta$roi <- as.integer(meta$roi)

#* Filter to target regions and ROIs
meta <- meta |>
  dplyr::filter(region %in% c("w32", "w214", "w271", "roiAD"),
                roi %in% c(1, 4, 8))

#* Loop through files
for (i in seq_len(nrow(meta))) {
  row <- meta[i, ]
  
  message("[INFO] Evaluating SPM for subject ", row$subjectID,
          ", region ", row$region, ", ROI ", row$roi, "...")
  
  # Define paths
  niftiFile     <- file.path(spmDir, paste0("binary_swww", row$subjectID, "_", row$region, "_", row$roi, ".nii"))
  roiMaskRegion <- ifelse(row$region == "roiAD", "wroiAD", row$region)
  roiFile       <- file.path(roiDir, paste0("wwwx", roiMaskRegion, "_redim_crop_squ_flipLR_newDim_", row$subjectID, ".nii"))
  outputDir     <- file.path(resultsDir, paste0("ROI", row$roi))
  outputCSV     <- file.path(outputDir, paste0("sens_esp_SPM_", row$region, "_", row$roi, ".csv"))
  
  if (!file.exists(niftiFile)) {
    warning("[WARNING] SPM file missing: ", niftiFile)
    next
  }
  if (!file.exists(roiFile)) {
    warning("[WARNING] ROI file missing: ", roiFile)
    next
  }
  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
  
  # Load points
  detectedPoints <- neuroSCC::getSPMbinary(niftiFile, paramZ = paramZ)
  truePoints <- neuroSCC::processROIs(
    roiFile = roiFile,
    region = roiMaskRegion,
    number = row$subjectID,
    save = FALSE,
    verbose = FALSE
  )
  
  # Skip if no true signal at z-slice
  sliceTruePoints <- subset(truePoints, z == paramZ & pet == 1)
  if (nrow(sliceTruePoints) == 0) {
    message("[WARNING] No ROI voxels at z = ", paramZ, " for ", row$subjectID, " — skipping.")
    next
  }
  
  # Evaluate
  subjectMetrics <- neuroSCC::calculateMetrics(
    detectedPoints = detectedPoints,
    truePoints     = truePoints,
    totalCoords    = dims,
    regionName     = paste0(row$region, "_", row$roi)
  )
  subjectMetrics$subject <- row$subjectID
  
  # Append to or create CSV
  if (!file.exists(outputCSV)) {
    readr::write_csv(subjectMetrics, outputCSV)
  } else {
    existing <- readr::read_csv(outputCSV, show_col_types = FALSE)
    combined <- dplyr::bind_rows(existing, subjectMetrics)
    readr::write_csv(combined, outputCSV)
  }
  
  rm(detectedPoints, truePoints, sliceTruePoints, subjectMetrics)
}

#* Clean up
rm(spmDir, resultsDir, roiDir, maskPath, dims, spmFiles, meta,
   outputDir, outputCSV, roiFile, niftiFile, roiMaskRegion)

### ==================================================== ###
###     7) SUMMARY TABLES (1 vs Group)                  ###
### ==================================================== ###

library(dplyr)
library(readr)
library(tidyr)

#* Paths
resultsDir <- file.path("Results", "z35", "results")
outputFileComplete <- file.path(resultsDir, "SCC_vs_SPM_complete.RDS")
outputFileSummary  <- file.path(resultsDir, "SCC_vs_SPM.RDS")

#* Initialize full dataset
SCC_vs_SPM_complete <- data.frame()

#* Define valid regions and ROIs to include
target_regions <- c("w32", "w214", "w271", "roiAD")
target_rois <- c(1, 4, 8)

#* Loop through region × ROI folders
for (region in target_regions) {
  for (roi in target_rois) {
    
    resultFolder <- file.path(resultsDir, paste0("ROI", roi))
    
    fileSCC <- file.path(resultFolder, paste0("sens_esp_SCC_", region, "_", roi, ".csv"))
    fileSPM <- file.path(resultFolder, paste0("sens_esp_SPM_", region, "_", roi, ".csv"))
    
    #* Process SCC file
    if (file.exists(fileSCC)) {
      tempSCC <- read_csv(fileSCC, show_col_types = FALSE)
      tempSCC <- tidyr::separate(tempSCC, col = region, into = c("region", "roi"), sep = "_", remove = FALSE)
      tempSCC <- mutate(tempSCC, method = "SCC", roi = as.integer(roi))
      tempSCC <- filter(tempSCC, region %in% target_regions, roi %in% target_rois)
      valid_cols <- c("method", "region", "roi", "sensitivity", "specificity", "PPV", "NPV", "subject")
      present_cols <- intersect(valid_cols, colnames(tempSCC))
      tempSCC <- tempSCC[, present_cols]
      
      SCC_vs_SPM_complete <- bind_rows(SCC_vs_SPM_complete, tempSCC)
    }
    
    #* Process SPM file
    if (file.exists(fileSPM)) {
      tempSPM <- read_csv(fileSPM, show_col_types = FALSE)
      tempSPM <- tidyr::separate(tempSPM, col = region, into = c("region", "roi"), sep = "_", remove = FALSE)
      tempSPM <- mutate(tempSPM, method = "SPM", roi = as.integer(roi))
      tempSPM <- filter(tempSPM, region %in% target_regions, roi %in% target_rois)
      valid_cols <- c("method", "region", "roi", "sensitivity", "specificity", "PPV", "NPV", "subject")
      present_cols <- intersect(valid_cols, colnames(tempSPM))
      tempSPM <- tempSPM[, present_cols]
      
      SCC_vs_SPM_complete <- bind_rows(SCC_vs_SPM_complete, tempSPM)
    }
  }
}

#* Compute summary (mean ± SD)
SCC_vs_SPM <- SCC_vs_SPM_complete |>
  group_by(method, region, roi) |>
  summarise(
    sensMEAN = mean(sensitivity, na.rm = TRUE),
    sensSD   = sd(sensitivity, na.rm = TRUE),
    espMEAN  = mean(specificity, na.rm = TRUE),
    espSD    = sd(specificity, na.rm = TRUE),
    ppvMEAN  = mean(PPV, na.rm = TRUE),
    ppvSD    = sd(PPV, na.rm = TRUE),
    npvMEAN  = mean(NPV, na.rm = TRUE),
    npvSD    = sd(NPV, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  arrange(method, region, roi)

#* Update region and ROI factor labels
region_levels <- c("w32", "w214", "w271", "roiAD")
region_labels <- c("ROI 1", "ROI 2", "ROI 3", "ROI 4")
roi_levels <- c(1, 4, 8)
roi_labels <- c("10", "40", "80")  # Hypoactivity levels (%)

SCC_vs_SPM_complete <- SCC_vs_SPM_complete %>%
  filter(region %in% region_levels, roi %in% roi_levels) %>%
  mutate(
    region = factor(region, levels = region_levels, labels = region_labels),
    roi = factor(roi, levels = roi_levels, labels = roi_labels)
  )

SCC_vs_SPM <- SCC_vs_SPM %>%
  filter(region %in% region_levels, roi %in% roi_levels) %>%
  mutate(
    region = factor(region, levels = region_levels, labels = region_labels),
    roi = factor(roi, levels = roi_levels, labels = roi_labels)
  )

#* Save summary table
saveRDS(SCC_vs_SPM, outputFileSummary)
write_csv(SCC_vs_SPM, gsub("\\.RDS$", ".csv", outputFileSummary))

#* Save full table
saveRDS(SCC_vs_SPM_complete, outputFileComplete)
write_csv(SCC_vs_SPM_complete, gsub("\\.RDS$", ".csv", outputFileComplete))

#* Clean up
rm(tempSCC, tempSPM, fileSCC, fileSPM, resultFolder,
   outputFileComplete, outputFileSummary, SCC_vs_SPM_complete, SCC_vs_SPM)

### ==================================================== ###
### 7B) STATISTICAL SIGNIFICANCE TESTS (Paired T-Tests) ###
### ==================================================== ###

library(dplyr)
library(tidyr)
library(readr)

#* Helper to assign stars
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("")
}

#* Ensure tibble structure
SCC_vs_SPM_complete <- as_tibble(SCC_vs_SPM_complete)

#* Initialize empty list
pvalue_table <- list()

#* Loop through region × roi
for (reg in levels(SCC_vs_SPM_complete$region)) {
  for (r in levels(SCC_vs_SPM_complete$roi)) {
    
    subset_data <- SCC_vs_SPM_complete %>%
      filter(region == reg, roi == r)
    
    # Skip if either method is missing
    if (!all(c("SCC", "SPM") %in% unique(subset_data$method))) next
    
    # Convert to wide format
    wide_data <- subset_data %>%
      dplyr::select(method, subject, sensitivity, specificity, PPV, NPV) %>%
      pivot_wider(names_from = method, values_from = c(sensitivity, specificity, PPV, NPV)) %>%
      drop_na()
    
    # Paired t-tests
    t_sens <- t.test(wide_data$sensitivity_SCC, wide_data$sensitivity_SPM, paired = TRUE)
    t_esp  <- t.test(wide_data$specificity_SCC, wide_data$specificity_SPM, paired = TRUE)
    t_ppv  <- t.test(wide_data$PPV_SCC, wide_data$PPV_SPM, paired = TRUE)
    t_npv  <- t.test(wide_data$NPV_SCC, wide_data$NPV_SPM, paired = TRUE)
    
    # Store results
    pvalue_table[[paste(reg, r, sep = "_")]] <- tibble(
      region = reg,
      roi = r,
      p_sens = t_sens$p.value,
      sig_sens = get_stars(t_sens$p.value),
      p_esp  = t_esp$p.value,
      sig_esp  = get_stars(t_esp$p.value),
      p_ppv  = t_ppv$p.value,
      sig_ppv = get_stars(t_ppv$p.value),
      p_npv  = t_npv$p.value,
      sig_npv = get_stars(t_npv$p.value)
    )
  }
}

#* Combine into single tibble
pvalue_table_1vsGroup <- bind_rows(pvalue_table)

#* Save outputs
write_csv(pvalue_table_1vsGroup, "Results/z35/results/pvalue_table_1vsGroup.csv")
saveRDS(pvalue_table_1vsGroup, "Results/z35/results/pvalue_table_1vsGroup.RDS")

#* Print preview
print(pvalue_table_1vsGroup)

#* Clean up environment
rm(get_stars, wide_data, subset_data, t_sens, t_esp, t_ppv, t_npv)


### ==================================================== ###
### 8) VISUALIZATIONS (1 vs Group)                      ###
### ==================================================== ###

#* Load evaluation and summary data
SCC_vs_SPM_complete <- readRDS("Results/z35/results/SCC_vs_SPM_complete.RDS")
SCC_vs_SPM <- readRDS("Results/z35/results/SCC_vs_SPM.RDS")

#* Load required packages
library(tidyverse)
library(lemon)
library(gridExtra)
library(ggridges)
library(viridis)
library(dotwhisker)

#* Set output folder for figures
setwd("Results/z35/Figures")

### ----------------------------------------------------------------
### Sensitivity, Specificity, PPV, NPV across all filtered regions with asterisks
### ----------------------------------------------------------------

plot_roi_grouped_with_asterisks <- function(table, pval_table, p_col = "p_sens", metric = "sensitivity", y_label = "Sensitivity (%)") {
  library(ggplot2)
  library(dplyr)
  
  dodge_width <- 0.8
  
  # 1. Prepare the input table
  table <- table %>%
    mutate(
      method = factor(method, levels = c("SCC", "SPM")),
      roi = factor(roi, levels = c("10", "40", "80")),
      region = factor(region)
    )
  
  # 2. Clean and label brackets (only asterisks)
  pval_clean <- pval_table %>%
    dplyr::select(region, roi, !!sym(p_col)) %>%
    mutate(
      region = as.character(region),
      roi = as.character(roi),
      group1 = "SCC",
      group2 = "SPM",
      label = case_when(
        !!sym(p_col) <= 0.001 ~ "***",
        !!sym(p_col) <= 0.01  ~ "**",
        !!sym(p_col) <= 0.05  ~ "*",
        TRUE                  ~ "ns"
      )
    )
  
  # 3. Match panel structure
  valid_panels <- table %>%
    distinct(region, roi) %>%
    mutate(across(everything(), as.character))
  pval_clean <- semi_join(pval_clean, valid_panels, by = c("region", "roi"))
  
  # 4. Compute height of bracket
  bracket_y <- table %>%
    group_by(region, roi) %>%
    summarise(y.position = max(.data[[metric]], na.rm = TRUE) + 5, .groups = "drop") %>%
    mutate(across(c(region, roi), as.character))
  
  # 5. Merge brackets
  bracket_data <- pval_clean %>%
    left_join(bracket_y, by = c("region", "roi")) %>%
    mutate(
      x = as.numeric(factor(roi, levels = c("10", "40", "80"))),
      x1 = x - dodge_width / 4,
      x2 = x + dodge_width / 4,
      label_y = y.position + 2
    )
  
  # 6. Final plot
  ggplot(table, aes(x = roi, y = .data[[metric]], fill = method)) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = 0.6,
      outlier.size = 0.5,
      lwd = 0.25
    ) +
    facet_wrap(~region, ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(
      limits = c(0, NA),
      breaks = c(25, 50, 75, 100),
      expand = expansion(mult = c(0, 0.1))
    ) +
    xlab("Hypoactivity (%)") + ylab(y_label) +
    theme_minimal(base_family = "serif") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 15),
      panel.spacing = unit(1, "lines")
    ) +
    # Bracket bar
    geom_segment(
      data = bracket_data,
      aes(x = x1, xend = x2, y = y.position, yend = y.position),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # Left vertical arm
    geom_segment(
      data = bracket_data,
      aes(x = x1, xend = x1, y = y.position, yend = y.position - 2),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # Right vertical arm
    geom_segment(
      data = bracket_data,
      aes(x = x2, xend = x2, y = y.position, yend = y.position - 2),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # Asterisk-only label
    geom_text(
      data = bracket_data,
      aes(x = x, y = label_y, label = label),
      size = 5,
      family = "serif",
      inherit.aes = FALSE
    )
}

# Set wd() for figure export
setwd("~/GitHub/PhD-2024-SCC-vs-SPM-SinglePatient-vs-Group/Results/z35/Figures")

# Create all 4 plots
graph_sens <- plot_roi_grouped_with_asterisks(
  table = SCC_vs_SPM_complete,
  pval_table = pvalue_table_1vsGroup,
  p_col = "p_sens",
  metric = "sensitivity",
  y_label = "Sensitivity (%)"
)

graph_esp <- plot_roi_grouped_with_asterisks(
  table = SCC_vs_SPM_complete,
  pval_table = pvalue_table_1vsGroup,
  p_col = "p_esp",
  metric = "specificity",
  y_label = "Specificity (%)"
)

graph_ppv <- plot_roi_grouped_with_asterisks(
  table = SCC_vs_SPM_complete,
  pval_table = pvalue_table_1vsGroup,
  p_col = "p_ppv",
  metric = "PPV",
  y_label = "Positive Predictive Value (%)"
)

graph_npv <- plot_roi_grouped_with_asterisks(
  table = SCC_vs_SPM_complete,
  pval_table = pvalue_table_1vsGroup,
  p_col = "p_npv",
  metric = "NPV",
  y_label = "Negative Predictive Value (%)"
)

# Save to PNGs
ggsave("sens_FILTERED.png", plot = graph_sens, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("esp_FILTERED.png", plot = graph_esp, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("ppv_FILTERED.png", plot = graph_ppv, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("npv_FILTERED.png", plot = graph_npv, width = 28.95, height = 18.3, units = "cm", dpi = 600)

### ----------------------------------------------------------------
### Sensitivity, Specificity, PPV, NPV across all filtered regions
### ----------------------------------------------------------------

plot_metric_by_region <- function(metric_name, y_label) {
  ggplot(data = SCC_vs_SPM_complete, aes(x = roi, y = .data[[metric_name]])) +
    geom_boxplot(aes(fill = method), outlier.size = 0.5, lwd = 0.25) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 100)) +
    xlab("Hypoactivity (%)") +
    ylab(y_label) +
    facet_wrap(~region, ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal(base_family = "serif") +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          legend.text = element_text(size = 14))
}

graph_sens <- plot_metric_by_region("sensitivity", "Sensitivity (%)")
graph_esp  <- plot_metric_by_region("specificity", "Specificity (%)")
graph_ppv  <- plot_metric_by_region("PPV", "Positive Predictive Value (%)")
graph_npv  <- plot_metric_by_region("NPV", "Negative Predictive Value (%)")

ggsave("sens_1vsGroup.png", plot = graph_sens, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("esp_1vsGroup.png",  plot = graph_esp,  width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("ppv_1vsGroup.png",  plot = graph_ppv,  width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("npv_1vsGroup.png",  plot = graph_npv,  width = 28.95, height = 18.3, units = "cm", dpi = 600)

combined_filtered <- grid.arrange(graph_sens, graph_esp, graph_ppv, graph_npv, ncol = 2)
ggsave("summary_metrics_1vsGroup.png", plot = combined_filtered, width = 38, height = 28, units = "cm", dpi = 600)

### ==================================================== ###
### OPTIONAL VISUALIZATION EXPERIMENTS (Ridge + Heatmap) ###
### ==================================================== ###

### -----------------------------------------------
### Ridge plot: Sensitivity by Region
### -----------------------------------------------

ridge_plot_sens <- ggplot(SCC_vs_SPM_complete, aes(x = sensitivity, y = region, fill = method)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_brewer(palette = "Set1") +
  theme_ridges(font_family = "serif") +
  labs(x = "Sensitivity (%)", y = "Region", title = "Density of Sensitivity by Region")

# ggsave("ridge_sensitivity_1vsGroup.png", ridge_plot_sens, width = 28, height = 18, units = "cm", dpi = 600)

### -----------------------------------------------
### Ridge plot: PPV by Region
### -----------------------------------------------

ridge_plot_ppv <- ggplot(SCC_vs_SPM_complete, aes(x = PPV, y = region, fill = method)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_brewer(palette = "Set1") +
  theme_ridges(font_family = "serif") +
  labs(x = "Positive Predictive Value (%)", y = "Region", title = "Density of PPV by Region")

# ggsave("ridge_ppv_1vsGroup.png", ridge_plot_ppv, width = 28, height = 18, units = "cm", dpi = 600)

### ----------------------------------------------------------
### Double-faceted Heatmap: Sensitivity & Specificity (1vsGroup)
### ----------------------------------------------------------

library(ggplot2)
library(patchwork)
library(grid)

# Heatmap for Sensitivity (1vsGroup)
heatmap_sens_facet <- ggplot(SCC_vs_SPM, aes(x = factor(roi), y = region, fill = sensMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = "Region",
    title = "Mean Sensitivity by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

# Heatmap for Specificity (1vsGroup)
heatmap_spec_facet <- ggplot(SCC_vs_SPM, aes(x = factor(roi), y = region, fill = espMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = NULL,
    title = "Mean Specificity by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm"),
    axis.title.y = element_blank()
  )

# Combine both with shared legend
combined_heatmap_sens_esp_1vsGroup <- heatmap_sens_facet + heatmap_spec_facet +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Show combined plot
combined_heatmap_sens_esp_1vsGroup

# Save to file
ggsave("combined_heatmap_sens_esp_1vsGroup.png", combined_heatmap_sens_esp_1vsGroup, width = 28, height = 20, units = "cm", dpi = 600)

### ----------------------------------------------------------
### Double-faceted Heatmap: PPV & NPV (1vsGroup)
### ----------------------------------------------------------

# Heatmap for PPV (1vsGroup)
heatmap_ppv_facet <- ggplot(SCC_vs_SPM, aes(x = factor(roi), y = region, fill = ppvMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = "Region",
    title = "Mean PPV by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

# Heatmap for NPV (1vsGroup)
heatmap_npv_facet <- ggplot(SCC_vs_SPM, aes(x = factor(roi), y = region, fill = npvMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = NULL,
    title = "Mean NPV by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm"),
    axis.title.y = element_blank()
  )

# Combine both with shared guide
combined_heatmap_ppv_npv_1vsGroup <- heatmap_ppv_facet + heatmap_npv_facet +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Show plot
combined_heatmap_ppv_npv_1vsGroup

# Save to file
ggsave("combined_heatmap_ppv_npv_1vsGroup.png", combined_heatmap_ppv_npv_1vsGroup, width = 28, height = 20, units = "cm", dpi = 600)


### -----------------------------------------------
### Double-faceted Heatmap: PPV
### -----------------------------------------------

heatmap_ppv_facet <- ggplot(SCC_vs_SPM, aes(x = factor(roi), y = region, fill = ppvMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_viridis(name = "PPV", option = "C", limits = c(0, 100)) +
  labs(x = "Hypoactivity Level (%)", y = "Region", title = "PPV (mean) by Method") +
  theme_minimal(base_family = "serif")

# ggsave("heatmap_ppv_1vsGroup.png", heatmap_ppv_facet, width = 28, height = 14, units = "cm", dpi = 600)



