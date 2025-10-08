# MGWNBR Analysis with Boruta Feature Selection
# Load required libraries
setwd("E:\\GWR SIMILARITY\\paper\\code\\mgwnbr\\R")
library(parallel)
library(MASS)
library(sp)
library(Boruta)
library(randomForest)

# Load MGWNBR function
source("E:\\GWR SIMILARITY\\paper\\code\\mgwnbr\\R\\mgwnbr4.R")

cat("=== MGWNBR Analysis with Boruta Feature Selection ===\n")
cat("Loading data_test1.csv...\n")

# Load dataset
data <- read.csv("E:\\GWR SIMILARITY\\paper\\code\\mgwnbr\\R\\data_test1.csv")
cat("Data loaded:", nrow(data), "observations,", ncol(data), "variables\n")

# Define variables
col_names <- names(data)
x_coord_col <- col_names[1]  # longitude column
y_coord_col <- col_names[2]  # latitude column
y_var <- col_names[3]        # dependent variable column
x_vars <- col_names[4:length(col_names)]  # independent variables

cat("Coordinates:", x_coord_col, "and", y_coord_col, "\n")
cat("Response variable:", y_var, "\n")
cat("Predictor variables:", length(x_vars), "variables\n")

# === STEP 1: BORUTA FEATURE SELECTION ===
cat("\n=== STEP 1: Running Boruta Feature Selection ===\n")

# Prepare data for Boruta
kde_vars <- x_vars  # All KDE variables
boruta_data <- data[, c(y_var, kde_vars)]

# Check for missing values
missing_count <- sum(is.na(boruta_data))
if (missing_count > 0) {
  cat("Warning: Found", missing_count, "missing values\n")
  boruta_data <- na.omit(boruta_data)
  cat("After removing missing values:", nrow(boruta_data), "observations\n")
}

# Run Boruta feature selection
cat("Running Boruta feature selection...\n")
set.seed(123)
boruta_result <- Boruta(crash ~ ., data = boruta_data, maxRuns = 200, doTrace = 2)

# Get selected features
confirmed_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
cat("\nBoruta Results:\n")
cat("Confirmed features:", length(confirmed_features), "\n")
cat("Selected features:", paste(confirmed_features, collapse=", "), "\n")

# Display Boruta summary
print(boruta_result)

# Save Boruta results
boruta_summary <- data.frame(
  Variable = confirmed_features,
  Decision = "Confirmed",
  MeanImp = apply(boruta_result$ImpHistory[, confirmed_features, drop = FALSE], 2, mean),
  MedianImp = apply(boruta_result$ImpHistory[, confirmed_features, drop = FALSE], 2, median)
)
write.csv(boruta_summary, "boruta_results.csv", row.names = FALSE)
cat("Boruta results saved as: boruta_results.csv\n")

# === STEP 2: MGWNBR WITH SELECTED FEATURES ===
cat("\n=== STEP 2: Running MGWNBR with Selected Features ===\n")

# Create filtered dataset with selected features
if (length(confirmed_features) > 0) {
  selected_data <- data[, c(x_coord_col, y_coord_col, y_var, confirmed_features)]
  cat("Using", length(confirmed_features), "selected features for MGWNBR\n")
} else {
  cat("Warning: No features selected by Boruta! Using all features.\n")
  selected_data <- data
  confirmed_features <- x_vars
}

# Build formula with selected features
formula_str <- paste(y_var, "~", paste(confirmed_features, collapse = " + "))
formula_obj <- as.formula(formula_str)
cat("Formula:", formula_str, "\n")

# Run MGWNBR with adaptive_bsq_smr method
cat("\nRunning MGWNBR with adaptive_bsq_smr method...\n")
start_time <- Sys.time()

result_adaptive_bsq_T <- mgwnbr(
  data = selected_data,
  formula = formula_obj,
  lat = y_coord_col,
  long = x_coord_col,
  method = "adaptive_bsq_smr",
  model = "negbin",
  mgwr = TRUE,
  bandwidth = "cv",
  globalmin = FALSE,
  verbose = TRUE
)

end_time <- Sys.time()
runtime <- end_time - start_time

cat("\n=== MGWNBR Analysis Completed ===\n")
cat("Runtime:", round(as.numeric(runtime), 2), "seconds\n")
cat("Selected features used:", length(confirmed_features), "\n")
cat("Features:", paste(confirmed_features, collapse=", "), "\n")

# Display key results
if (!is.null(result_adaptive_bsq_T)) {
  cat("\n=== FINAL RESULTS ===\n")
  if (!is.null(result_adaptive_bsq_T$measures)) {
    cat("Final AIC:", round(result_adaptive_bsq_T$measures$AIC, 2), "\n")
    cat("Final AICc:", round(result_adaptive_bsq_T$measures$AICc, 2), "\n")
    cat("Final Deviance:", round(result_adaptive_bsq_T$measures$deviance, 2), "\n")
    if (!is.null(result_adaptive_bsq_T$measures$percent_deviance)) {
      cat("Percent Deviance Explained:", round(result_adaptive_bsq_T$measures$percent_deviance, 2), "%\n")
    }
  }
  if (!is.null(result_adaptive_bsq_T$band)) {
    cat("Final Bandwidths:", paste(round(result_adaptive_bsq_T$band, 2), collapse=", "), "\n")
  }
  if (!is.null(result_adaptive_bsq_T$alphaW_values)) {
    cat("Final AlphaW Values:", paste(round(result_adaptive_bsq_T$alphaW_values, 3), collapse=", "), "\n")
  }
  
  # Display detailed iteration results if available
  if (!is.null(result_adaptive_bsq_T$iteration_results)) {
    cat("\n=== ITERATION PERFORMANCE SUMMARY ===\n")
    iteration_results <- result_adaptive_bsq_T$iteration_results
    if (is.data.frame(iteration_results)) {
      print(iteration_results)
    } else {
      cat("Iteration results available but not in expected format\n")
    }
  }
  
  # Display parameter estimates summary
  if (!is.null(result_adaptive_bsq_T$mgwr_param_estimates)) {
    cat("\n=== PARAMETER ESTIMATES SUMMARY ===\n")
    param_estimates <- result_adaptive_bsq_T$mgwr_param_estimates
    if (is.data.frame(param_estimates)) {
      cat("Parameter estimates summary:\n")
      print(summary(param_estimates))
    }
  }
}

cat("\nAnalysis completed successfully!\n")

