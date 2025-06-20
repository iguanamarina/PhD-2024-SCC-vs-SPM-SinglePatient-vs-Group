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
### 7B) STATISTICAL SIGNIFICANCE TESTS with testCompareR ###
### ==================================================== ###

# NOTE for future users:
# This section uses compareR() to assess statistical significance in a 1vsGroup setting.
# It compares SCC vs SPM for Sensitivity, Specificity, PPV, and NPV.
# File structures:
#   SCC: Results/z35/1vsGroup/SCC/SCC_C{n}_{region}_{roi}.RData
#   SPM: Results/z35/1vsGroup/SPM/binary_swwwC{n}_{region}_{roi}.nii

library(testCompareR)
library(dplyr)
library(readr)

# Load custom triplet builder
source("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/Contrastes de Hipótesis/generateCompareRTriplets.R")

# Define paths
base_path <- "Results/z35/1vsGroup"
scc_dir <- file.path(base_path, "SCC")
spm_dir <- file.path(base_path, "SPM")
mask_path <- "Auxiliary Files/new_mask.nii"

# Load grid
dims <- neuroSCC::getDimensions(mask_path)
grid <- expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c("x", "y")]

# Regions and ROIs (filtered)
regions <- c("w32", "w214", "w271", "roiAD")  # order matters
rois <- c(1, 4, 8)

# Get list of subject IDs from SCC filenames
scc_files <- list.files(scc_dir, pattern = "^SCC_C\\d+_.*\\.RData$", full.names = TRUE)
subjects <- sort(unique(as.numeric(sub(".*SCC_C(\\d+)_.*", "\\1", scc_files))))

# Storage
pvalResults <- list()

# Loop region × roi
for (region in regions) {
  for (roi in rois) {
    
    message("[INFO] Region ", region, " | ROI ", roi)
    tripletList <- list()
    
    for (subject in subjects) {
      subj_str <- paste0("C", subject)
      
      # Paths
      scc_file <- file.path(scc_dir, paste0("SCC_", subj_str, "_", region, "_", roi, ".RData"))
      spm_file <- file.path(spm_dir, paste0("binary_swww", subj_str, "_", region, "_", roi, ".nii"))
      
      if (!file.exists(scc_file) || !file.exists(spm_file)) {
        warning("[WARNING] Missing files for subject ", subj_str, " | region ", region, " | ROI ", roi)
        next
      }
      
      env <- new.env()
      load(scc_file, envir = env)
      scc_obj <- env$SCC_1vsG
      scc_points <- neuroSCC::getPoints(scc_obj)$positivePoints
      spm_points <- neuroSCC::getSPMbinary(spm_file, paramZ = paramZ)
      
      # ROI mask
      roi_region <- ifelse(region == "roiAD", "wroiAD", region)
      roi_file <- file.path("roisNormalizadas", paste0("wwwx", roi_region, "_redim_crop_squ_flipLR_newDim_", subj_str, ".nii"))
      if (!file.exists(roi_file)) {
        warning("[WARNING] Missing ROI file for subject ", subj_str)
        next
      }
      
      true_points <- neuroSCC::processROIs(roi_file, region, subj_str, save = FALSE, verbose = FALSE)
      true_slice <- subset(true_points, z == paramZ & pet == 1, select = c("x", "y"))
      if (nrow(true_slice) == 0) next
      
      triplet <- generateCompareRTriplets(
        grid = grid,
        sccCoords = scc_points,
        spmCoords = spm_points,
        roiCoords = true_slice
      )
      tripletList[[length(tripletList) + 1]] <- triplet
    }
    
    if (length(tripletList) == 0) {
      warning("[WARNING] No valid subjects for ", region, " ROI ", roi)
      next
    }
    
    allTriplets <- bind_rows(tripletList) %>%
      mutate(across(everything(), as.numeric))
    
    # Run test
    result <- compareR(
      df = allTriplets,
      test1 = "SCC_pred",
      test2 = "SPM_pred",
      gold = "ROI_truth",
      interpret = FALSE,
      multi_corr = "holm",
      alpha = 0.05,
      sesp = TRUE,
      ppvnpv = TRUE,
      plrnlr = TRUE,
      test.names = c("SCC", "SPM"),
      dp = 2
    )
    
    # Extract
    row <- tibble(
      region = region,
      roi = roi,
      sens_SCC = result$acc$accuracies$SCC["Sensitivity", "Estimate"],
      se_sens_SCC = result$acc$accuracies$SCC["Sensitivity", "SE"],
      sens_SPM = result$acc$accuracies$SPM["Sensitivity", "Estimate"],
      se_sens_SPM = result$acc$accuracies$SPM["Sensitivity", "SE"],
      p_sens = result$acc$sens.p.adj,
      
      spec_SCC = result$acc$accuracies$SCC["Specificity", "Estimate"],
      se_spec_SCC = result$acc$accuracies$SCC["Specificity", "SE"],
      spec_SPM = result$acc$accuracies$SPM["Specificity", "Estimate"],
      se_spec_SPM = result$acc$accuracies$SPM["Specificity", "SE"],
      p_spec = result$acc$spec.p.adj,
      
      ppv_SCC = result$pv$predictive.values$SCC["PPV", "Estimate"],
      se_ppv_SCC = result$pv$predictive.values$SCC["PPV", "SE"],
      ppv_SPM = result$pv$predictive.values$SPM["PPV", "Estimate"],
      se_ppv_SPM = result$pv$predictive.values$SPM["PPV", "SE"],
      p_ppv = result$pv$ppv.p.adj,
      
      npv_SCC = result$pv$predictive.values$SCC["NPV", "Estimate"],
      se_npv_SCC = result$pv$predictive.values$SCC["NPV", "SE"],
      npv_SPM = result$pv$predictive.values$SPM["NPV", "Estimate"],
      se_npv_SPM = result$pv$predictive.values$SPM["NPV", "SE"],
      p_npv = result$pv$npv.p.adj,
      
      lrpos_SCC = result$lr$likelihood.ratios$SCC["PLR", "Estimate"],
      se_lrpos_SCC = result$lr$likelihood.ratios$SCC["PLR", "SE"],
      lrpos_SPM = result$lr$likelihood.ratios$SPM["PLR", "Estimate"],
      se_lrpos_SPM = result$lr$likelihood.ratios$SPM["PLR", "SE"],
      p_lrpos = result$lr$plr.p.adj,
      
      lrneg_SCC = result$lr$likelihood.ratios$SCC["NLR", "Estimate"],
      se_lrneg_SCC = result$lr$likelihood.ratios$SCC["NLR", "SE"],
      lrneg_SPM = result$lr$likelihood.ratios$SPM["NLR", "Estimate"],
      se_lrneg_SPM = result$lr$likelihood.ratios$SPM["NLR", "SE"],
      p_lrneg = result$lr$nlr.p.adj
    )
    
    pvalResults[[length(pvalResults) + 1]] <- row
  }
}

# Output
pvalTable_1vsGroup <- bind_rows(pvalResults)

# Mapear orden correcto de regiones
region_order <- c("w32", "w214", "w271", "roiAD")
pvalTable_1vsGroup <- pvalTable_1vsGroup %>%
  filter(region %in% region_order, roi %in% c(1, 4, 8)) %>%
  mutate(region = factor(region, levels = region_order))

# Significance stars
getStars <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("")
}

pvalTable_1vsGroup <- pvalTable_1vsGroup %>%
  mutate(
    sig_sens = sapply(p_sens, getStars),
    sig_spec = sapply(p_spec, getStars),
    sig_ppv  = sapply(p_ppv, getStars),
    sig_npv  = sapply(p_npv, getStars),
    sig_lrpos = sapply(p_lrpos, getStars),
    sig_lrneg = sapply(p_lrneg, getStars)
  )

# Reorder columns for output
pvalTable_1vsGroup <- pvalTable_1vsGroup %>%
  dplyr::select(
    region, roi,
    sens_SCC, se_sens_SCC, sens_SPM, se_sens_SPM, p_sens, sig_sens,
    spec_SCC, se_spec_SCC, spec_SPM, se_spec_SPM, p_spec, sig_spec,
    ppv_SCC, se_ppv_SCC, ppv_SPM, se_ppv_SPM, p_ppv, sig_ppv,
    npv_SCC, se_npv_SCC, npv_SPM, se_npv_SPM, p_npv, sig_npv,
    lrpos_SCC, se_lrpos_SCC, lrpos_SPM, se_lrpos_SPM, p_lrpos, sig_lrpos,
    lrneg_SCC, se_lrneg_SCC, lrneg_SPM, se_lrneg_SPM, p_lrneg, sig_lrneg
  )

# Save
write_csv(pvalTable_1vsGroup, file = file.path(base_path, "pvalue_table_compareR.csv"))
saveRDS(pvalTable_1vsGroup, file = file.path(base_path, "pvalue_table_compareR.RDS"))
# pvalue_table_compareR <- readRDS("~/GitHub/PhD-2024-SCC-vs-SPM-SinglePatient-vs-Group/Results/z35/1vsGroup/pvalue_table_compareR.RDS")


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
### Sensitivity, Specificity, PPV, NPV — Single vs Group — with Asterisks
### ----------------------------------------------------------------

plot_metric_boxplot <- function(pval_data,
                                metric = c("sensitivity", "specificity", "ppv", "npv"),
                                reps = 25) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  metric <- match.arg(metric)
  dodge_width <- 0.8
  
  # -----------------------------
  # Configuración por métrica
  # -----------------------------
  metric_map <- list(
    sensitivity = "sens",
    specificity = "spec",
    ppv         = "ppv",
    npv         = "npv"
  )
  raw_name <- metric_map[[metric]]
  value_col_SCC <- paste0(raw_name, "_SCC")
  value_col_SPM <- paste0(raw_name, "_SPM")
  se_col_SCC    <- paste0("se_", raw_name, "_SCC")
  se_col_SPM    <- paste0("se_", raw_name, "_SPM")
  p_col         <- paste0("p_", raw_name)
  
  y_label_name <- switch(metric,
                         sensitivity = "Sensitivity (%)",
                         specificity = "Specificity (%)",
                         ppv         = "Positive Predictive Value (%)",
                         npv         = "Negative Predictive Value (%)"
  )
  
  if (metric == "ppv") {
    y_limits <- c(0, 40)
    y_breaks <- seq(0, 40, 5)
    bracket_offset <- 4
    rabito <- 1
    texto_offset <- 1
  } else if (metric == "sensitivity") {
    y_limits <- c(0, 108)
    y_breaks <- seq(0, 100, by = 20)
    bracket_offset <- 6
    rabito <- 2
    texto_offset <- 1
  } else if (metric == "npv") {
    y_limits <- c(80, 110)  # << se sube ligeramente el techo
    y_breaks <- seq(80, 100, by = 20)
    bracket_offset <- 4
    rabito <- 1.5
    texto_offset <- 1
  } else {
    y_limits <- c(0, 100)
    y_breaks <- seq(0, 100, by = 20)
    bracket_offset <- 4
    rabito <- 1.5
    texto_offset <- 1
  }
  
  # -----------------------------
  # Expandir simulaciones
  # -----------------------------
  long_df <- pval_data %>%
    mutate(
      region = factor(region, levels = c("w32", "w214", "w271", "roiAD"),
                      labels = c("ROI 1", "ROI 2", "ROI 3", "ROI 4")),
      roi = factor(roi, levels = c(1, 4, 8), labels = c("10", "40", "80"))
    ) %>%
    rowwise() %>%
    mutate(
      value_SCC = list({
        sd_val <- .data[[se_col_SCC]]
        sd_val <- if (is.na(sd_val) || !is.finite(sd_val)) 1e-6 else max(sd_val, 1e-6)
        rnorm(reps, mean = .data[[value_col_SCC]], sd = sd_val)
      }),
      value_SPM = list({
        sd_val <- .data[[se_col_SPM]]
        sd_val <- if (is.na(sd_val) || !is.finite(sd_val)) 1e-6 else max(sd_val, 1e-6)
        rnorm(reps, mean = .data[[value_col_SPM]], sd = sd_val)
      })
    ) %>%
    unnest(cols = c(value_SCC, value_SPM), names_sep = "_") %>%
    pivot_longer(cols = c(value_SCC, value_SPM),
                 names_to = "method", values_to = "value",
                 names_pattern = "value_(.*)") %>%
    ungroup() %>%
    filter(is.finite(value))
  
  if (nrow(long_df) == 0) {
    stop(paste("No valid data for plotting", metric, "- all filtered out."))
  }
  
  valid_panels <- long_df %>%
    count(region, roi) %>%
    filter(n > 0) %>%
    dplyr::select(-n)
  
  # Etiquetas de significancia (NA → ns)
  pval_clean <- pval_data %>%
    mutate(
      region = factor(region, levels = c("w32", "w214", "w271", "roiAD"),
                      labels = c("ROI 1", "ROI 2", "ROI 3", "ROI 4")),
      roi = factor(roi, levels = c(1, 4, 8), labels = c("10", "40", "80")),
      label = case_when(
        is.na(!!rlang::sym(p_col))         ~ "ns",
        !!rlang::sym(p_col) <= 0.001 ~ "***",
        !!rlang::sym(p_col) <= 0.01  ~ "**",
        !!rlang::sym(p_col) <= 0.05  ~ "*",
        TRUE                         ~ "ns"
      )
    ) %>%
    dplyr::select(region, roi, label)
  
  # Altura de brackets
  bracket_y <- long_df %>%
    group_by(region, roi) %>%
    summarise(ymax = max(value, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(ymax)) %>%
    mutate(
      y.position = pmin(ymax + bracket_offset, y_limits[2] - 1),
      y.bottom   = y.position - rabito,
      label_y    = pmin(y.position + texto_offset + 1, y_limits[2])  # 👈 extra espacio para "ns"
    )
  
  bracket_data <- pval_clean %>%
    semi_join(bracket_y, by = c("region", "roi")) %>%
    left_join(bracket_y, by = c("region", "roi")) %>%
    semi_join(valid_panels, by = c("region", "roi")) %>%
    mutate(
      x = as.numeric(roi),
      x1 = x - dodge_width / 4,
      x2 = x + dodge_width / 4
    )
  
  long_df <- semi_join(long_df, valid_panels, by = c("region", "roi"))
  
  # -----------------------------
  # Plot final
  # -----------------------------
  ggplot(long_df, aes(x = roi, y = value, fill = method)) +
    geom_boxplot(
      aes(color = method),
      position = position_dodge(width = dodge_width),
      width = 0.6,
      outlier.size = 0.5,
      lwd = 0.25
    ) +
    facet_wrap(~region, ncol = 2) +
    scale_fill_manual(values = c("SCC" = "#d73027", "SPM" = "#4575b4")) +
    scale_color_manual(values = c("SCC" = "#d73027", "SPM" = "#4575b4")) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Hypoactivity (%)", y = y_label_name) +
    theme_minimal(base_family = "serif") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 15),
      panel.spacing = unit(1, "lines"),
      legend.title = element_blank()
    ) +
    geom_segment(data = bracket_data, aes(x = x1, xend = x2, y = y.position, yend = y.position),
                 linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(data = bracket_data, aes(x = x1, xend = x1, y = y.position, yend = y.bottom),
                 linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(data = bracket_data, aes(x = x2, xend = x2, y = y.position, yend = y.bottom),
                 linewidth = 0.5, inherit.aes = FALSE) +
    geom_text(data = bracket_data, aes(x = x, y = label_y, label = label),
              size = 5, family = "serif", inherit.aes = FALSE)
}



# --------------------------------------------
# Crear y guardar gráficos para Single vs Group
# --------------------------------------------

#* Set output folder for figures
setwd("Results/z35/Figures")

# Create Graphics
graph_sens_1vs <- plot_metric_boxplot(pvalue_table_compareR, metric = "sensitivity")
graph_spec_1vs <- plot_metric_boxplot(pvalue_table_compareR, metric = "specificity")
graph_ppv_1vs  <- plot_metric_boxplot(pvalue_table_compareR, metric = "ppv")
graph_npv_1vs  <- plot_metric_boxplot(pvalue_table_compareR, metric = "npv")

# Save as PNGs
png("sens_1vsGroup_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_sens_1vs)
dev.off()

png("esp_1vsGroup_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_spec_1vs)
dev.off()

png("ppv_1vsGroup_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_ppv_1vs)
dev.off()

png("npv_1vsGroup_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_npv_1vs)
dev.off()



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

### ----------------------------------------------------------
### Preparar datos: SCC_vs_SPM para 1vsGroup
### ----------------------------------------------------------

SCC_vs_SPM <- pvalue_table_compareR %>%
  mutate(
    region = factor(region, levels = c("w32", "w214", "w271", "roiAD"),
                    labels = c("ROI 1", "ROI 2", "ROI 3", "ROI 4")),
    roi = factor(roi, levels = c(1, 4, 8), labels = c("10", "40", "80"))
  ) %>%
  dplyr::select(region, roi,
                sens_SCC, sens_SPM,
                spec_SCC, spec_SPM,
                ppv_SCC, ppv_SPM,
                npv_SCC, npv_SPM) %>%
  pivot_longer(cols = starts_with("sens_"),
               names_to = "method_sens", names_prefix = "sens_",
               values_to = "sensMEAN") %>%
  pivot_longer(cols = starts_with("spec_"),
               names_to = "method_spec", names_prefix = "spec_",
               values_to = "espMEAN") %>%
  pivot_longer(cols = starts_with("ppv_"),
               names_to = "method_ppv", names_prefix = "ppv_",
               values_to = "ppvMEAN") %>%
  pivot_longer(cols = starts_with("npv_"),
               names_to = "method_npv", names_prefix = "npv_",
               values_to = "npvMEAN") %>%
  filter(method_sens == method_spec,
         method_sens == method_ppv,
         method_sens == method_npv) %>%
  rename(method = method_sens) %>%
  dplyr::select(region, roi, method, sensMEAN, espMEAN, ppvMEAN, npvMEAN)


### ----------------------------------------------------------
### Double-faceted Heatmap: Sensitivity & Specificity (1vsGroup)
### ----------------------------------------------------------

library(ggplot2)
library(patchwork)
library(grid)

heatmap_sens_facet <- ggplot(SCC_vs_SPM, aes(x = roi, y = region, fill = sensMEAN)) +
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
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

heatmap_spec_facet <- ggplot(SCC_vs_SPM, aes(x = roi, y = region, fill = espMEAN)) +
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
    y = "",
    title = "Mean Specificity by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

combined_heatmap_sens_esp_1vsGroup <- heatmap_sens_facet + heatmap_spec_facet +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# Mostrar y guardar
combined_heatmap_sens_esp_1vsGroup
ggsave("combined_heatmap_sens_esp_1vsGroup.png", combined_heatmap_sens_esp_1vsGroup,
       width = 28, height = 20, units = "cm", dpi = 600)

### ----------------------------------------------------------
### Double-faceted Heatmap: PPV & NPV (1vsGroup)
### ----------------------------------------------------------

heatmap_ppv_facet <- ggplot(SCC_vs_SPM, aes(x = roi, y = region, fill = ppvMEAN)) +
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
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

heatmap_npv_facet <- ggplot(SCC_vs_SPM, aes(x = roi, y = region, fill = npvMEAN)) +
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
    y = "",
    title = "Mean NPV by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

combined_heatmap_ppv_npv_1vsGroup <- heatmap_ppv_facet + heatmap_npv_facet +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# Mostrar y guardar
combined_heatmap_ppv_npv_1vsGroup
ggsave("combined_heatmap_ppv_npv_1vsGroup.png", combined_heatmap_ppv_npv_1vsGroup,
       width = 28, height = 20, units = "cm", dpi = 600)


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



