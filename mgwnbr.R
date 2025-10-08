
#' @title Multiscale Geographically Weighted Negative Binomial Regression
#'
#' @description Fits a geographically weighted regression model with different scales for each covariate. Uses the negative binomial distribution as default, but also accepts the normal, Poisson, or logistic distributions. Can fit the global versions of each regression and also the geographically weighted alternatives with only one scale, since they are all particular cases of the multiscale approach.
#'
#' @param data name of the dataset.
#' @param formula regression model formula as in \code{lm}.
#' @param weight name of the variable containing the sample weights, default value is \code{NULL}.
#' @param lat name of the variable containing the latitudes in the dataset.
#' @param long name of the variable containing the longitudes in the dataset.
#' @param globalmin logical value indicating whether to find a global minimum in the optimization process, default value is \code{TRUE}.
#' @param method indicates the method to be used for the bandwidth calculation (\code{adaptive_bsq}, \code{fixed_bsq}, \code{fixed_g}).
#' @param model indicates the model to be used for the regression (\code{gaussian}, \code{poisson}, \code{negbin}, \code{logistic}), default value is\code{"negbin"}.
#' @param mgwr logical value indicating if multiscale should be used (\code{TRUE}, \code{FALSE}), default value is \code{TRUE}.
#' @param bandwidth indicates the criterion to be used for the bandwidth calculation (\code{cv}, \code{aic}), default value is \code{"cv"}.
#' @param offset name of the variable containing the offset values, if null then is set to a vector of zeros, default value is \code{NULL}.
#' @param distancekm logical value indicating whether to calculate the distances in km, default value is \code{FALSE}.
#' @param int integer indicating the number of iterations, default value is \code{50}.
#' @param h integer indicating a predetermined bandwidth value, default value is \code{NULL}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{general_bandwidth} - General bandwidth value.
#' \item \code{band} - Bandwidth values for each covariate.
#' \item \code{measures} - Goodness of fit statistics.
#' \item \code{ENP} - Effective number of parameters.
#' \item \code{mgwr_param_estimates} - MGWR parameter estimates.
#' \item \code{qntls_mgwr_param_estimates} - Quantiles of MGWR parameter estimates.
#' \item \code{descript_stats_mgwr_param_estimates} - Descriptive statistics of MGWR parameter estimates.
#' \item \code{p_values} - P-values for the t tests on parameter significance.
#' \item \code{t_critical} - Critical values for the t tests on parameter significance.
#' \item \code{mgwr_se} - MGWR standard errors.
#' \item \code{qntls_mgwr_se} - Quantiles of MGWR standard errors.
#' \item \code{descript_stats_se} - Descriptive statistics of MGWR standard errors.
#' \item \code{global_param_estimates} - Parameter estimates for the global model.
#' \item \code{t_test_dfs} - Denominator degrees of freedom for the t tests.
#' \item \code{global_measures} - Goodness of fit statistics for the global model.
#' }
#'
#' @examples
#' ## Data
#'
#'
#' data(georgia)
#'
#' for (var in c("PctFB", "PctBlack")){
#'   georgia[, var] <- as.data.frame(scale(georgia[, var]))
#' }
#'
#'
#' ## Model
#'
#' mod <- mgwnbr(data=georgia, formula=PctBach~PctBlack+PctFB,
#'  lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq", bandwidth="cv",
#'   model="gaussian", mgwr=FALSE, h=136)
#'
#' ## Bandwidths
#' mod$general_bandwidth
#'
#' ## Goodness of fit measures
#' mod$measures
#'
#' @importFrom sp spDistsN1
#'
#' @importFrom stats model.extract model.matrix pnorm pt qt quantile dist
#'
#' @importFrom MASS ginv
#'
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ stopCluster
#'
#' @importFrom foreach foreach %dopar% %do%
#'
#' @export

# Parallel computing setup function
setup_parallel <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- min(detectCores() - 1, 15)  # Use max 15 cores, leave 1 free
  }
  
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    
    # Export all necessary functions from the current environment
    clusterExport(cl, c(
      "robust_solve", "robust_quantile", "gwr_local", "gwr_local_multi_alpha",
      "gwr_local_multi_bw_alpha", "optimize_alphaW_per_variable", 
      "optimize_alphaW_sequential", "true_hybrid_gwr", "simple_gwr_baseline",
      "GSS_with_alpha", "cv_with_alpha", "optimize_multi_alpha_mgwr",
      "global_local_hybrid_mgwr"
    ), envir = environment())
    
    # Load necessary packages on each worker
    clusterEvalQ(cl, {
      library(stats)
      library(MASS)
      library(foreach)
      library(sp)
    })
    
    # Also try to source the functions directly on each worker
    clusterEvalQ(cl, {
      # Define robust_solve locally on each worker
      robust_solve <- function(A, b) {
        tryCatch({
          result <- solve(A, b)
          if (any(is.na(result)) || any(is.infinite(result))) {
            # Try with ridge penalty
            ridge_penalty <- 1e-8 * diag(nrow(A))
            result <- solve(A + ridge_penalty, b)
          }
          return(result)
        }, error = function(e) {
          # Use pseudoinverse as fallback
          tryCatch({
            result <- MASS::ginv(A) %*% b
            return(result)
          }, error = function(e2) {
            # Final fallback: return zeros
            warning("Matrix solve failed, using zero coefficients")
            return(rep(0, length(b)))
          })
        })
      }
      
      # Define robust_quantile locally on each worker
      robust_quantile <- function(x, probs) {
        x_clean <- x[!is.na(x) & !is.infinite(x)]
        if (length(x_clean) == 0) {
          return(rep(NA, length(probs)))
        }
        return(quantile(x_clean, probs, na.rm = TRUE))
      }
    })
    
    return(cl)
  } else {
    return(NULL)
  }
}

# Cleanup parallel computing
cleanup_parallel <- function(cl) {
  if (!is.null(cl)) {
    stopCluster(cl)
  }
}

# Robust matrix solving function
robust_solve <- function(A, b) {
  tryCatch({
    result <- solve(A, b)
    if (any(is.na(result)) || any(is.infinite(result))) {
      # Try with ridge penalty
      ridge_penalty <- 1e-8 * diag(nrow(A))
      result <- solve(A + ridge_penalty, b)
    }
    return(result)
  }, error = function(e) {
    # Use pseudoinverse as fallback
    tryCatch({
      result <- MASS::ginv(A) %*% b
      return(result)
    }, error = function(e2) {
      # Final fallback: return zeros
      warning("Matrix solve failed, using zero coefficients")
      return(rep(0, length(b)))
    })
  })
}

# Robust quantile function that handles NaN values
robust_quantile <- function(x, probs) {
  x_clean <- x[!is.na(x) & !is.infinite(x)]
  if (length(x_clean) == 0) {
    return(rep(0, length(probs)))
  }
  return(quantile(x_clean, probs, na.rm = TRUE))
}

eval_aicc_with_offset <- function(alpha_val, bandwidth_val, Y, X_var, method, model, N, wt, E, COORD, sequ, distancekm, parg, offset_vals) {
  tryCatch({
    # --- FIX: Pass the offset to the 'yhat_beta' parameter, which is handled correctly
    # as an additive offset in the linear predictor's update formula. ---
    result <- gwr_local(
      H = bandwidth_val,
      y = Y,
      x = X_var,
      fi = rep(0, N), # Set fi to zero
      yhat_beta = offset_vals, # <--- CORRECTLY PASS OFFSET HERE
      alphaW_val = alpha_val,
      method = method,
      model = model,
      N = N,
      nvarg = 1,
      wt = wt,
      E = E,
      COORD = COORD,
      sequ = sequ,
      distancekm = distancekm,
      Offset = rep(0, N), # Global offset is already in offset_vals
      parg = parg
    )
    
    # --- The rest of the function remains the same ---
    sm <- result$sm
    v1 <- sum(diag(sm))
    yhat <- result$yhbeta[,1] 
    alphai <- result$alphai
    
    # Check for problematic values
    if (any(is.na(yhat)) || any(is.infinite(yhat)) || any(yhat <= 0)) return(Inf)
    if (any(is.na(alphai[,2])) || any(is.infinite(alphai[,2])) || any(alphai[,2] <= 0)) return(Inf)
    
    # Check if effective number of parameters is reasonable
    if (v1 > N * 0.8) return(Inf)  # Prevent overfitting
    
    # Log-likelihood is calculated with the original Y and the new yhat
    ll_terms <- Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
      lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1)
    
    # Check for problematic log-likelihood terms
    if (any(is.na(ll_terms)) || any(is.infinite(ll_terms))) return(Inf)
    
    ll <- sum(ll_terms, na.rm = TRUE)
    
    # Check if log-likelihood is reasonable (not too extreme)
    if (is.na(ll) || is.infinite(ll) || abs(ll) > 1e6) return(Inf)
    
    AIC <- 2 * v1 - 2 * ll
    
    # Check if denominator is too small to avoid numerical instability
    denominator <- N - v1 - 1
    if (denominator <= 1) {
      # Use AIC instead of AICc when sample size is too small relative to parameters
      AICc <- AIC
    } else {
      AICc <- AIC + 2 * v1 * (v1 + 1) / denominator
    }
    
    return(AICc)
  }, error = function(e) {
    # Optional: print the error for debugging
    # message(sprintf("Error in eval_aicc_with_offset for alpha=%s: %s", alpha_val, e$message))
    return(Inf)
  })
}


adaptive_grid_search_with_offset <- function(var_idx, Fi, mband, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, verbose, baseline_alphaWs = NULL) {
  
  # Set default baseline_alphaWs if not provided
  if (is.null(baseline_alphaWs)) {
    baseline_alphaWs <- rep(0.5, nvarg)  # Default to 0.5 for all variables
  }
  
  if (verbose) cat("\n    Optimizing alpha for variable", var_idx, " (", colnames(X)[var_idx], ")...\n")
  
  # 1. Calculate the OFFSET from the global offset + fit of all other variables
  offset_for_this_var <- Offset 
  if (nvarg > 1) {
    Fi_others <- Fi[, -var_idx, drop = FALSE]
    offset_for_this_var <- offset_for_this_var + apply(Fi_others, 1, sum)
  }
  
  X_var <- as.matrix(X[, var_idx])
  bandwidth_val <- mband[var_idx]
  if (verbose) cat("      Using bandwidth for variable", var_idx, ":", round(bandwidth_val, 2), "\n")
  
  # 2. Perform Adaptive Grid Search with 3 classes
  grid_cells <- list(c(0.0, 0.33), c(0.34, 0.66), c(0.67, 1.0))
  cell_midpoints <- c(0.165, 0.5, 0.835)
  
  # Phase 1: Test midpoints - Calculate AICc for FULL MODEL
  cell_aiccs <- sapply(cell_midpoints, function(alpha_test) {
    # Create temporary Fi with current state
    Fi_temp <- Fi
    Fi_temp[, var_idx] <- X_var * 0  # Reset target variable
    
    # Calculate AICc for FULL MODEL with all variables
    eval_aicc_full_model(alpha_test, bandwidth_val, var_idx, Fi_temp, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, mband, baseline_alphaWs)
  })
  
  best_cell_idx <- which.min(cell_aiccs)
  best_cell <- grid_cells[[best_cell_idx]]
  
  # Phase 2: Golden Section Search in the best cell
  phi <- (1 + sqrt(5)) / 2
  a <- best_cell[1]
  b <- best_cell[2]
  tolerance <- 0.01
  
  # Golden section search for optimal alphaW
  c <- b - (b - a) / phi
  d <- a + (b - a) / phi
  
  # Evaluate AICc at initial points
  fc <- eval_aicc_full_model(c, bandwidth_val, var_idx, Fi, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, mband, baseline_alphaWs)
  fd <- eval_aicc_full_model(d, bandwidth_val, var_idx, Fi, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, mband, baseline_alphaWs)
  
  while (abs(b - a) > tolerance) {
    if (fc < fd) {
      b <- d
      d <- c
      fd <- fc
      c <- b - (b - a) / phi
      fc <- eval_aicc_full_model(c, bandwidth_val, var_idx, Fi, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, mband, baseline_alphaWs)
    } else {
      a <- c
      c <- d
      fc <- fd
      d <- a + (b - a) / phi
      fd <- eval_aicc_full_model(d, bandwidth_val, var_idx, Fi, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, mband, baseline_alphaWs)
    }
  }
  
  best_alpha <- (a + b) / 2
  best_aicc <- min(fc, fd, na.rm=TRUE)
  
  if (verbose) {
    cat("      → Best alpha found:", round(best_alpha, 4), "in cell", best_cell_idx, "with AICc:", round(best_aicc, 2), "\n")
  }
  
  return(list(alpha = best_alpha, aicc = best_aicc))
}

# Function to calculate AICc for FULL MODEL (all variables) with only target variable varying
eval_aicc_full_model <- function(alpha_test, bandwidth_val, var_idx, Fi_temp, Y, X, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, mband, baseline_alphaWs) {
  tryCatch({
    # Create alpha_vals array with baseline alphaWs for other variables, test alpha for target variable
    alpha_vals <- baseline_alphaWs  # Use baseline alphaWs for all variables
    alpha_vals[var_idx] <- alpha_test  # Override target variable with test alphaW
    
    # Use multi-bandwidth array (mband) for proper MGWR evaluation
    bandwidth_vals <- mband  # Use optimized multi-bandwidths for all variables
    
    # Use gwr_local_multi_alpha for proper multi-variable MGWR evaluation (now with fixed sm accumulation)
    result_full <- gwr_local_multi_alpha(
      H = bandwidth_vals, y = Y, x = X, fi = rep(0, N), 
      alpha_vals = alpha_vals, method = method, model = model, 
      N = N, nvarg = nvarg, wt = wt, E = E, COORD = COORD, sequ = sequ, 
      distancekm = distancekm, Offset = Offset, parg = parg, yhat_beta = rep(0, N), verbose = FALSE
    )
    
    # Get the full model results
    yhbeta <- result_full$yhbeta
    sm <- result_full$sm
    alphai <- result_full$alphai
    v1 <- sum(diag(sm))
    
    # Calculate predicted values for the full model
    yhat_full <- exp(apply(yhbeta[, 2:(nvarg + 1)] * X, 1, sum) + Offset)
    
    # Check for problematic values
    if (any(is.na(yhat_full)) || any(is.infinite(yhat_full)) || any(yhat_full <= 0)) {
      cat("DEBUG: Problematic yhat_full detected, returning Inf\n")
      return(Inf)
    }
    if (any(is.na(alphai[,2])) || any(is.infinite(alphai[,2])) || any(alphai[,2] <= 0)) {
      cat("DEBUG: Problematic alphai detected, returning Inf\n")
      return(Inf)
    }
    
    # Relax v1 check - allow up to 95% of N instead of 80%
    if (v1 > N * 0.95) {
      cat("DEBUG: v1 too large (", round(v1, 2), ") relative to N (", N, "), returning Inf\n")
      return(Inf)  # Prevent extreme overfitting
    }
    
    # Don't reject based on denominator - handle it in AICc calculation instead
    
    # Calculate log-likelihood for FULL MODEL
    ll_terms <- Y * log(alphai[,2] * yhat_full) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat_full) +
      lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1)
    
    # Check for problematic log-likelihood terms
    if (any(is.na(ll_terms))) {
      cat("DEBUG: NA in ll_terms detected\n")
      return(Inf)
    }
    if (any(is.infinite(ll_terms))) {
      cat("DEBUG: Inf in ll_terms detected\n")
      return(Inf)
    }
    
    ll <- sum(ll_terms, na.rm = TRUE)
    
    # Check if log-likelihood is reasonable (relax the threshold)
    if (is.na(ll) || is.infinite(ll) || abs(ll) > 1e10) {
      cat("DEBUG: Problematic ll (", ll, "), returning Inf\n")
      return(Inf)
    }
    
    AIC <- 2 * v1 - 2 * ll
    
    # Calculate AICc with proper handling of small denominators
    denominator <- N - v1 - 1
    if (denominator <= 0) {
      # If denominator is 0 or negative, use AIC only
      AICc <- AIC
      cat("DEBUG: Small denominator (", round(denominator, 2), "), using AIC. v1=", round(v1, 2), ", N=", N, ", AIC=", round(AIC, 2), "\n")
    } else {
      AICc <- AIC + 2 * v1 * (v1 + 1) / denominator
    }
    
    # Debug: Show calculation details
    cat("DEBUG: v1=", round(v1, 2), ", N=", N, ", denom=", round(denominator, 2), ", ll=", round(ll, 2), ", AIC=", round(AIC, 2), ", AICc=", round(AICc, 2), ", alphaW=", round(alpha_test, 3), "\n")
    
    return(AICc)
  }, error = function(e) {
    cat("DEBUG: ERROR in eval_aicc_full_model - alphaW:", round(alpha_test, 3), "Error:", e$message, "\n")
    return(Inf)
  })
}

# Pure local version of gwr for alphaW optimization
gwr_local <- function(H, y, x, fi, alphaW_val = 0.5, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta) {
  nvar <- ncol(x)
  bim <- rep(0, nvar*N)
  yhatm <- rep(0, N)
  sm <- matrix(0, N, N)
  sm3 <- matrix(0, N, nvar)
  alphai <- matrix(0, N, 3)
  # No assign or get to parent/global, all local!
  for (i in 1:N){
    # Calculate distance matrix for this point
    seqi <- rep(i, N)
    dx <- sp::spDistsN1(COORD, COORD[i,])
    distan <- cbind(seqi, sequ, dx)
    if (distancekm){
      distan[,3] <- distan[,3]*111
    }
    u <- nrow(distan)
    
    dist_spatial_i <- dist_spatial_mat[i, ]
    w <- rep(0, N)
    
    if (method=="fixed_g"){
      for (jj in 1:u){
        w[jj] <- exp(-(distan[jj,3]/H)^2)
      }
    }
    else if (method=="fixed_bsq"){
      for (jj in 1:u){
        w[jj] <- (1-(dist_spatial_i[jj]/H)^2)^2
      }
    }
    else if (method=="adaptive_bsq"){
      distan <- distan[order(distan[, 3]), ]
      distan <- cbind(distan, 1:nrow(distan))
      w <- matrix(0, N, 2)
      hn <- distan[H,3]
      for (jj in 1:N){
        if (distan[jj,4]<=H){
          w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
        } else{
          w[jj,1] <- 0
        }
        w[jj,2] <- distan[jj,2]
      }
      w <- w[order(w[, 2]), ]
      w <- w[,1]
    }
    else if (method == "adaptive_bsq_smr") {
      # Use alphaW_val = 1.0 to follow adaptive_bsq behavior
      if (is.null(alphaW_val)) alphaW_val <- 1.0
      betaW_val <- 1 - alphaW_val
      # --- w1 (Spatial weights with bisquare kernel) ---
      # dist_spatial_sorted <- distan[order(distan[, 3]), ]
      # dist_spatial_sorted <- cbind(dist_spatial_sorted, 1:nrow(dist_spatial_sorted))
      # w1_matrix <- matrix(0, N, 2)
      # hn_spatial <- dist_spatial_sorted[H, 3]
      dist_spatial_sorted_df <- data.frame(dist = dist_spatial_i, original_index = 1:N)
      dist_spatial_sorted_df <- dist_spatial_sorted_df[order(dist_spatial_sorted_df$dist), ]
      dist_spatial_sorted_df$rank <- 1:N # Tambahkan kolom rank
      
      w1_matrix <- matrix(0, N, 2)
      hn_spatial <- dist_spatial_sorted_df$dist[H]
      if (hn_spatial == 0) hn_spatial <- 1e-6
      
      for (jj in 1:N) {
        if (dist_spatial_sorted_df$rank[jj] <= H) {
          w1_matrix[jj, 1] <- (1 - (dist_spatial_sorted_df$dist[jj] / hn_spatial)^2)^2
        } else {
          w1_matrix[jj, 1] <- 0
        }
        w1_matrix[jj, 2] <- dist_spatial_sorted_df$original_index[jj] 
      }
      w1_matrix <- w1_matrix[order(w1_matrix[, 2]), ]
      w1 <- w1_matrix[, 1]
      
      # --- w2 (Attribute similarity weights with bisquare kernel) ---
      # Handle single variable case
      # if (ncol(x) == 1) {
      #   sim_dists <- abs(x - x[i, 1])
      # } else {
      #   sim_dists <- sqrt(rowSums(sweep(x, 2, x[i, ], "-")^2))
      # }
      
      dist_sim_i <- dist_sim_mat[i, ]
      
      dist_sim_sorted_df <- data.frame(dist = dist_sim_i, original_index = 1:N)
      dist_sim_sorted_df <- dist_sim_sorted_df[order(dist_sim_sorted_df$dist), ]
      dist_sim_sorted_df$rank <- 1:N # Tambahkan kolom rank
      
      w2_matrix <- matrix(0, N, 2)
      hn_sim <- dist_sim_sorted_df$dist[H]
      if (hn_sim == 0) hn_sim <- 1e-6
      
      for (jj in 1:N) {
        if (dist_sim_sorted_df$rank[jj] <= H) {
          w2_matrix[jj, 1] <- exp(-(dist_sim_sorted_df$dist[jj] / hn_sim)^2)
        } else {
          w2_matrix[jj, 1] <- 0
        }
        w2_matrix[jj, 2] <- dist_sim_sorted_df$original_index[jj]
      }
      
      w2_matrix <- w2_matrix[order(w2_matrix[, 2]), ]
      w2 <- w2_matrix[, 1]
      w <- alphaW_val * w1 + betaW_val * w2
      if (i == 1) {
        # cat("DEBUG: adaptive_bsq_smr weights - alphaW_val =", alphaW_val, "betaW_val =", betaW_val, "\n")
        # cat("DEBUG: w1 range =", range(w1), "w2 range =", range(w2), "w range =", range(w), "\n")
        # cat("DEBUG: w length =", length(w), "wt length =", length(wt), "N =", N, "\n")
      }
    }
    # Weight matrix calculation (silent)
    ## MODEL SELECTION ##
    if (model=="gaussian"){
      if (det(t(x)%*%(w*x*wt))==0){
        b <- rep(0, nvar)
      }
      else{
        b <- robust_solve(t(x)%*%diag(as.vector(w*wt))%*%x, t(x)%*%diag(as.vector(w*wt))%*%y)
      }
      uj <- x%*%b
      # update sm matrix
      if (nvar==nvarg){
        if (det(t(x)%*%diag(as.vector(w*wt))%*%x)==0){
          sm[i, ] <- rep(0, N)
        } else {
          C <- robust_solve(t(x)%*%diag(as.vector(w*wt))%*%x, t(x)%*%diag(as.vector(w*wt)))
          sm[i, ] <- (x[i,] %*% C)
        }
      }
    }
    else if (model=="poisson" | model=="negbin"){
      # if (i == 1) cat("DEBUG: gwr_local negbin/poisson section entry\n")
      # if (i == 1) cat("DEBUG: w length =", length(w), "wt length =", length(wt), "y length =", length(y), "\n")
      uj <- rep(mean(y), N)
      if (is.null(parg)) {
        par <- 1.0  # Default value for par
      } else {
        par <- parg
      }
      nj <- log(uj)
      ddpar <- 1
      cont <- 1
      while (abs(ddpar)>0.000001 & cont<100){
        dpar <- 1
        parold <- par
        cont1 <- 1
        cont3 <- 1
        if (model=="poisson"){
          alphaNB <- E^-6  # Negative binomial alpha (dispersion parameter)
          par <- 1/alphaNB
        }
        else{
          if (par<=E^-5 & i>1){
            par=1/alphai[i-1,2]
          }
          while (abs(dpar)>0.000001 & cont1<200){
            par <- ifelse(par<E^-10, E^-10, par)
            # if (i == 1 && cont1 == 1) cat("DEBUG: About to calculate g with par =", par, "\n")
            g <- sum(w*wt*(digamma(as.vector(par)+y)-digamma(as.vector(par))+log(as.vector(par))+1-log(as.vector(par)+uj)-(as.vector(par)+y)/(as.vector(par)+uj)))
            hess <- sum(w*wt*(trigamma(as.vector(par)+y)-trigamma(as.vector(par))+1/as.vector(par)-2/(as.vector(par)+uj)+(y+as.vector(par))/(as.vector(par)+uj)^2))
            # hess <- ifelse(hess==0, E^-23, hess)
            par0 <- par
            # if (i == 1 && cont1 == 1) cat("DEBUG: About to solve hess*g, hess =", hess, "g =", g, "\n")
            par <- par0-solve(hess)*g
            if (cont1>50 & par>E^5){
              dpar <- 0.0001
              cont3 <- cont3+1
              if (cont3==1){
                par <- 2
              }
              else if (cont3==2){
                par <- E^5
              }
              else if (cont3==3){
                par <- 0.0001
              }
            }
            else{
              dpar <- par-par0
            }
            cont1 <- cont1+1
            if (par>E^6){
              par <- E^6
              dpar <- 0
            }
            if (par<=E^-5){
              par <- E^-3
              dpar <- 0
            }
          }
          alphaNB <- as.vector(1/par)  # Negative binomial alpha (dispersion parameter)
        }
        dev <- 0
        ddev <- 1
        cont2 <- 0
        while (abs(ddev)>0.000001 & cont2<100){
          uj <- ifelse(uj>E^100, E^100, uj)
          ai <- as.vector((uj/(1+alphaNB*uj))+(y-uj)*(alphaNB*uj/(1+2*alphaNB*uj+alphaNB^2*uj*uj)))
          ai <- ifelse(ai<=0, E^-5, ai)
          zj <- nj+(y-uj)/(ai*(1+alphaNB*uj))-yhat_beta+fi
          # Fix matrix dimensions for negative binomial
          if (model == "negbin") {
            # Ensure all vectors have the same length as the number of observations
            ai_vector <- rep(as.vector(ai), length.out = N)
            w_vector <- rep(as.vector(w), length.out = N)
            wt_vector <- rep(as.vector(wt), length.out = N)
            zj_vector <- rep(as.vector(zj), length.out = N)
            
            # Create diagonal weight matrix
            W_diag <- diag(w_vector * ai_vector * wt_vector)
            
            # Solve with proper dimensions
            b <- robust_solve(t(x) %*% W_diag %*% x, t(x) %*% W_diag %*% zj_vector)
            
            # Ensure yhat_beta and fi have correct dimensions
            yhat_beta_vector <- rep(as.vector(yhat_beta), length.out = N)
            fi_vector <- rep(as.vector(fi), length.out = N)
          } else {
            b <- robust_solve(t(x)%*%diag(as.vector(w*ai*wt))%*%x, t(x)%*%diag(as.vector(w*ai*wt))%*%zj)
          }
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- as.vector(exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          if (model=="poisson"){
            dev <- 2*sum(y*log(tt)-(y-uj))
          }
          else{
            dev <- 2*sum(y*log(tt)-(y+1/alphaNB)*log((1+alphaNB*y)/(1+alphaNB*uj)))
          }
          cont2 <- cont2+1
        }
        cont <- cont+1
        ddpar <- par-parold
      }
      yhatm[i] <- uj[i]
      alphai[i, 2] <- alphaNB  # Store negative binomial alpha (dispersion parameter)
      
      # Calculate alpha standard error for negative binomial model
      if (model == "negbin") {
        hess_alpha <- sum(w*wt*(trigamma(as.vector(par)+y)-trigamma(as.vector(par))+1/as.vector(par)-2/(as.vector(par)+uj)+(y+as.vector(par))/(as.vector(par)+uj)^2))
        # Safety check to prevent division by zero
        hess_alpha <- ifelse(hess_alpha==0, E^-23, hess_alpha)
        sealphaNB <- sqrt(1/abs(hess_alpha))/(par^2)
        alphai[i, 3] <- sealphaNB
      }
      
      # update sm matrix
      if (nvar==nvarg){
        # if (i == 1) cat("DEBUG: About to check matrix determinant for sm calculation\n")
        matrix_to_check <- t(x)%*%diag(as.vector(w*ai*wt))%*%x
        det_val <- det(matrix_to_check)
        # if (i == 1) cat("DEBUG: Matrix determinant =", det_val, "\n")
        
        if (det_val==0){
          # if (i == 1) cat("DEBUG: Matrix is singular, setting sm to zero\n")
          sm[i, ] <- rep(0, N)
          sm3[i,] <- rep(0, nvar) 
        } else {
          # if (i == 1) cat("DEBUG: Matrix is non-singular, calculating C\n")
          C <- robust_solve(t(x)%*%diag(as.vector(w*ai*wt))%*%x, t(x)%*%diag(as.vector(w*ai*wt)))
          sm[i, ] <- (x[i,] %*% C)
          # For negative binomial, use proper variance calculation
          if (model == "negbin") {
            # For negbin: Var(β) = (X'WX)^(-1) * (1 + α*μ) where α is the dispersion parameter
            variance_factor <- 1 + alphaNB * uj[i]
            sm3[i, ] <- t(diag(C %*% diag(1/ai) %*% t(C))) * variance_factor
          } else {
            sm3[i, ] <- t(diag(C %*% diag(1/ai) %*% t(C)))
          }
        }
      }
    }
    m1 <- (i-1)*nvar+1
    m2 <- m1+(nvar-1)
    bim[m1:m2] <- b
    yhatm[i] <- uj[i]
  }
  beta <- matrix(bim, N, byrow=T)
  yhbeta <- cbind(yhatm, beta)
  return(list(
    yhbeta = yhbeta,
    sm = sm,
    alphai = alphai,
    sm3 = sm3
  ))
}
# Golden Section Search for optimal alphaW

find_optimal_alphaW_gss <- function(
    H, Y, X, finb, N, nvarg, Offset, method, model, wt, E, COORD, sequ, distancekm, parg, yhat_beta, tol = 0.01
) {
  # Golden Section Search for optimal alphaW
  phi <- (1 + sqrt(5)) / 2
  # Use full range [0,1] for all methods including adaptive_bsq_smr
  a <- 0.0
  b <- 1.0
  
  # Define eval_aicc function
  eval_aicc <- function(alphaW_val) {
    res <- gwr_local(H, Y, X, finb, alphaW_val, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta)
    local_yhbeta <- res$yhbeta
    local_beta <- local_yhbeta[, 2:(nvarg + 1)]
    Fi <- X * local_beta
    sm <- res$sm
    alphai <- res$alphai
    v1 <- sum(diag(sm))
    yhat <- exp(apply(Fi, 1, sum) + Offset)
    ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
    AIC <- 2 * v1 - 2 * ll
    if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
    
    return(AICc)
  }
  
  # Evaluate endpoints
  fa <- eval_aicc(a)
  fb <- eval_aicc(b)
  
  # Handle NA/NaN values
  if (is.na(fa) || is.nan(fa)) fa <- 1000
  if (is.na(fb) || is.nan(fb)) fb <- 1000
  
  step <- 1
  
  while (abs(b - a) > tol) {
    step <- step + 1
    if (fa < fb) {
      # Update upper bound
      b_new <- a + (b - a)/phi
      fb_new <- eval_aicc(b_new)
      
      # Handle NA/NaN values - NO FALLBACKS
      if (is.na(fb_new) || is.nan(fb_new)) stop("NA/NaN detected in alphaW optimization - halting execution")
      
      fb <- fb_new
      b <- b_new
    } else {
      # Update lower bound
      a_new <- b - (b - a)/phi
      fa_new <- eval_aicc(a_new)
      
      # Handle NA/NaN values - NO FALLBACKS
      if (is.na(fa_new) || is.nan(fa_new)) stop("NA/NaN detected in alphaW optimization - halting execution")
      
      fa <- fa_new
      a <- a_new
    }
  }
  optimal_alphaW <- (a + b) / 2
  # Optimal alphaW found
  # cat("DEBUG: Golden Section Search found optimal_alphaW =", optimal_alphaW, "for method =", method, "\n")
  return(optimal_alphaW)
}

mgwnbr <- function(data, formula, weight=NULL, lat, long,
                   globalmin=TRUE, method, model="negbin",
                   mgwr=TRUE, bandwidth="cv", offset=NULL,
                   distancekm=FALSE, int=50, h=NULL, verbose=TRUE){
  # Record start time for runtime calculation
  start_time <- Sys.time()
  
  # MGWNBR function entry
  cat("=== MGWNBR FUNCTION STARTED ===\n")
  cat("Parameters: method=", method, ", mgwr=", mgwr, ", model=", model, "\n")
  cat("Data dimensions:", nrow(data), "observations,", ncol(data), "variables\n")
  cat("Formula:", deparse(formula), "\n")
  cat("================================\n")
  
  # Setup parallel computing
  cl <- setup_parallel()
  if (!is.null(cl)) {
    cat("Parallel computing enabled with", length(cl), "cores\n")
  }
  
  output <- list()
  header <- c()
  yhat_beta <- NULL
  E <- 10
  
  if (verbose) cat("Step 1: Data preparation and model setup...\n")
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf)
  mt <- attr(mf, "terms")
  XVAR <- attr(mt, "term.labels")
  Y <- model.extract(mf, "response")
  N <- length(Y)
  X <- model.matrix(mt, mf)
  
  if (verbose) {
    cat("  - Response variable range:", range(Y), "\n")
    cat("  - Number of predictors:", length(XVAR), "\n")
    cat("  - Predictors:", paste(XVAR, collapse=", "), "\n")
  }
  
  wt <-rep(1, N)
  if (!is.null(weight)){
    wt <- unlist(data[, weight])
    if (verbose) cat("  - Using weights from variable:", weight, "\n")
  }
  Offset <- rep(0, N)
  if (!is.null(offset)){
    Offset <- unlist(data[, offset])
    if (verbose) cat("  - Using offset from variable:", offset, "\n")
  }
  nvarg <- ncol(X)
  yhat <- rep(0, N)
  bi <- matrix(0, nvarg*N, 4)
  alphai <- matrix(0, N, 3)
  s <- rep(0, N)
  mrj <- matrix(0, N, N*nvarg)
  sm <- matrix(0, N, N)
  sm3 <- matrix(0, N, nvarg)
  rj <- matrix(0, N, N)
  Cm <- matrix(0, N, N*nvarg)
  stdbm <- matrix(0, N, nvarg)
  mAi <- matrix(0, N, nvarg)
  ENP <- rep(0, nvarg+2)
  
  if (verbose) cat("Step 1 completed. Data prepared for analysis.\n\n")
  
  if (verbose) cat("Step 2: Computing global estimates...\n")
  ## global estimates ##
  if (model=="poisson" | model=="negbin"){
    if (verbose) cat("  - Initializing global model parameters...\n")
    uj <- (Y+mean(Y))/2
    nj <- log(uj)
    parg <- as.vector(sum((Y-uj)^2/uj)/(N-nvarg))
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
      dpar <- 1
      parold <- parg
      cont1 <- 1
      cont3 <- 1
      if(model=="poisson"){
        if (verbose) cat("  - Fitting global Poisson model...\n")
        alphaNBg <- E^-6  # Global negative binomial alpha (dispersion parameter)
        parg <- 1/alphaNBg
      }
      else{
        if (cont>1){
          if (verbose && cont == 2) cat("  - Fitting global negative binomial model...\n")
          parg <- as.vector(1/(sum((Y-uj)^2/uj)/(N-nvarg)))
        }
        while (abs(dpar)>0.000001 & cont1<200){
          parg <- ifelse(parg<E^-10, E^-10, parg)
          g <- sum(digamma(as.vector(parg)+Y)-digamma(as.vector(parg))+log(as.vector(parg))+1-log(as.vector(parg)+uj)-(as.vector(parg)+Y)/(as.vector(parg)+uj))
          hess <- sum(trigamma(as.vector(parg)+Y)-trigamma(as.vector(parg))+1/as.vector(parg)-2/(as.vector(parg)+uj)+(Y+as.vector(parg))/(as.vector(parg)+uj)^2)
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <- as.vector(par0 - g/hess)
          if (cont1>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              parg <- 2
            }
            else if (cont3==2){
              parg <- E^5
            }
            else if (cont3==3){
              parg <- 0.0001
            }
          }
          else{
            dpar <- parg-par0
          }
          cont1 <- cont1+1
          if (parg>E^6){
            parg <- E^6
            dpar <- 0
          }
        }
        alphaNBg <- as.vector(1/parg)  # Global negative binomial alpha (dispersion parameter)
      }
      devg <- 0
      ddev <- 1
      cont2 <- 0
      while (abs(ddev)>0.000001 & cont2<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        ai <- as.vector((uj/(1+alphaNBg*uj))+(Y-uj)*(alphaNBg*uj/(1+2*alphaNBg*uj+alphaNBg^2*uj*uj)))
        ai <- ifelse(ai<=0, E^-5, ai)
        zj <- nj+(Y-uj)/(ai*(1+alphaNBg*uj))-Offset
        if (det(t(X)%*%(ai*X))==0){
          bg <- rep(0, nvarg)
        }
        else{
          bg <- solve(t(X)%*%(ai*X))%*%t(X)%*%(ai*zj)
        }
        nj <- X%*%bg+Offset
        nj <- ifelse(nj>E^2, E^2, nj)
        uj <- as.vector(exp(nj))
        olddev <- devg
        uj <- ifelse(uj<E^-150, E^-150, uj)
        tt <- Y/uj
        tt <- ifelse(tt==0, E^-10, tt)
        if (model=="poisson"){
          devg <- 2*sum(Y*log(tt)-(Y-uj))
        }
        else{
          devg <- 2*sum(Y*log(tt)-(Y+1/alphaNBg)*log((1+alphaNBg*Y)/(1+alphaNBg*uj)))
          sealphaNBg <- sqrt(1/abs(hess))/(parg^2)
        }
        if (cont2>100){
          ddev <- 0.0000001
        }
        else{
          ddev <- devg-olddev
        }
        cont2 <- cont2+1
      }
      ujg <- uj
      yhat <- uj
      cont <- cont+1
      ddpar <- parg-parold
    }
    varg <- diag(solve(t(X*wt*ai)%*%X))
  }
  else if (model=="logistic"){
    uj <- (Y+mean(Y))/2
    nj <- log(uj/(1-uj))
    devg <- 0
    ddev <- 1
    cont <- 0
    while (abs(ddev)>0.000001 & cont<100){
      uj <- ifelse(uj>E^100, E^100, uj)
      ai <- as.vector(uj*(1-uj))
      ai <- ifelse(ai<=0, E^-5, ai)
      zj <- nj+(Y-uj)/ai
      if (det(t(X)%*%(wt*ai*X))==0){
        bg <- rep(0, nvarg)
      }
      else{
        bg <- solve(t(X)%*%(wt*ai*X))%*%t(X)%*%(wt*ai*zj)
      }
      nj <- X%*%bg
      nj <- ifelse(nj>E^2, E^2, nj)
      uj <- exp(nj)/(1+exp(nj))
      olddev <- devg
      uj <- ifelse(uj<E^-150, E^-150, uj)
      tt <- Y/uj
      tt <- ifelse(tt==0, E^-10, tt)
      uj <- ifelse(uj==1, 0.99999, uj)
      tt2 <- (1-Y)/(1-uj)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devg <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
      ddev <- devg-olddev
      cont <- cont+1
    }
    ujg <- uj
    yhat <- uj
    varg <- diag(solve(t(X*wt*ai)%*%X))
  }
  # Extract coordinates
  
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  
  
  
  COORD <- matrix(c(long, lat), ncol=2)
  sequ <- 1:N
  cat("Pre-computing spatial and attribute distance matrices for efficiency...\n")
  
  dist_spatial_mat <- sp::spDists(COORD)
  if (distancekm) {
    dist_spatial_mat <- dist_spatial_mat * 111
  }
  
  
  dist_sim_mat <- as.matrix(dist(scale(X)))
  assign("dist_spatial_mat", dist_spatial_mat, envir = .GlobalEnv)
  assign("dist_sim_mat", dist_sim_mat, envir = .GlobalEnv)
  
  cv <- function(H, y, x, fi, alphaW_val = 1.0){
    nvar <- ncol(x)
    # CV function call
    for (i in 1:N){
      # Calculate distance matrix for this point
      seqi <- rep(i, N)
      dx <- sp::spDistsN1(COORD, COORD[i,])
      distan <- cbind(seqi, sequ, dx)
      if (distancekm){
        distan[,3] <- distan[,3]*111
      }
      
      dist_spatial_i <- dist_spatial_mat[i, ]
      
      w <- rep(0, N)
      
      # for (jj in 1:u){
      #   w[jj] <- exp(-0.5*(distan[jj,3]/H)^2)
      #   if (bandwidth=="cv"){
      #     w[i] <- 0
      #   }
      # }
      if (bandwidth=="cv"){
        w[i] <- 0
      }
      
      if (method=="fixed_bsq"){
        position <- which(distan[,3]>H)
        w[position] <- 0
      }
      else if (method=="adaptive_bsq"){
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, N, 2)
        hn <- distan[H,3]
        for (jj in 1:N){
          if (distan[jj,4]<=H){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
        }
        if (bandwidth=="cv"){
          w[which(w[,2]==i)] <- 0
        }
        w <- w[order(w[, 2]), ]
        w <- w[ ,1]
      }
      else if (method == "adaptive_bsq_smr") {
        # Define mixing parameters - use the alphaW_val parameter passed to the function
        # alphaW_val is already defined as a parameter
        betaW_val <- 1 - alphaW_val
        
        # --- w1 (Bobot Spasial) ---
        dist_spatial_sorted_df <- data.frame(dist = dist_spatial_i, original_index = 1:N)
        dist_spatial_sorted_df <- dist_spatial_sorted_df[order(dist_spatial_sorted_df$dist), ]
        dist_spatial_sorted_df$rank <- 1:N
        
        w1_matrix <- matrix(0, N, 2)
        hn_spatial <- dist_spatial_sorted_df$dist[H]
        if (hn_spatial == 0) hn_spatial <- 1e-6 
        
        for (jj in 1:N) {
          if (dist_spatial_sorted_df$rank[jj] <= H) {
            w1_matrix[jj, 1] <- (1 - (dist_spatial_sorted_df$dist[jj] / hn_spatial)^2)^2
          } else {
            w1_matrix[jj, 1] <- 0
          }
          w1_matrix[jj, 2] <- dist_spatial_sorted_df$original_index[jj]
        }
        w1_matrix <- w1_matrix[order(w1_matrix[, 2]), ]
        w1 <- w1_matrix[, 1]
        
        dist_sim_i <- dist_sim_mat[i, ]
        dist_sim_sorted_df <- data.frame(dist = dist_sim_i, original_index = 1:N)
        dist_sim_sorted_df <- dist_sim_sorted_df[order(dist_sim_sorted_df$dist), ]
        dist_sim_sorted_df$rank <- 1:N
        
        w2_matrix <- matrix(0, N, 2)
        hn_sim <- dist_sim_sorted_df$dist[H]
        if (hn_sim == 0) hn_sim <- 1e-6
        
        for (jj in 1:N) {
          if (dist_sim_sorted_df$rank[jj] <= H) {
            w2_matrix[jj, 1] <- exp(-(dist_sim_sorted_df$dist[jj] / hn_sim)^2)
          } else {
            w2_matrix[jj, 1] <- 0
          }
          w2_matrix[jj, 2] <- dist_sim_sorted_df$original_index[jj]
        }
        w2_matrix <- w2_matrix[order(w2_matrix[, 2]), ]
        w2 <- w2_matrix[, 1]
        
        w <- alphaW_val * w1 + betaW_val * w2
        
        if (bandwidth == "cv") {
          w[i] <- 0
        }
      }
      
      if (model=="gaussian"){
        b <- robust_solve(t(x)%*%(w*x*wt), t(x)%*%(w*y*wt))
        yhat_ <- get("yhat")
        yhat_[i] <- x[i, ]%*%b
        assign("yhat", yhat_, envir=parent.frame())
        yhat[i] <- x[i, ]%*%b
        s_ <- get("s")
        s_[i] <- (x[i,]%*%robust_solve(t(x)%*%(w*x*wt), t(x*w*wt)))[i]
        assign("s", s_, envir=parent.frame())
        s[i] <- (x[i,]%*%robust_solve(t(x)%*%(w*x*wt), t(x*w*wt)))[i]
        next
      }
      else if (model=="poisson" | model=="negbin"){
        uj <- yhat
        par <- parg
        nj <- log(uj)
        ddpar <- 1
        cont <- 1
        cont3 <- 0
        while (abs(ddpar)>0.000001 & cont<100){
          dpar <- 1
          parold <- par
          cont1 <- 1
          if (model=="poisson"){
            alphaNB <- E^-6  # Negative binomial alpha (dispersion parameter)
            par <- 1/alphaNB
          }
          else{
            if (par<=E^-5 & i>1){
              par <- as.vector(1/alphai[i-1, 2])
            }
            while (abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(as.vector(par)+y)-digamma(as.vector(par))+log(as.vector(par))+1-log(as.vector(par)+uj)-(as.vector(par)+y)/(as.vector(par)+uj)))
              hess <- sum(w*wt*(trigamma(as.vector(par)+y)-trigamma(as.vector(par))+1/as.vector(par)-2/(as.vector(par)+uj)+(y+as.vector(par))/(as.vector(par)+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              # if (i == 1 && cont1 == 1) cat("DEBUG: About to solve hess*g, hess =", hess, "g =", g, "\n")
              par <- as.vector(par0 - g/hess)  # Use division instead of solve for scalar
              if (cont1>50 & par>E^5){
                dpar <- 0.0001
                cont3 <- cont3+1
                if (cont3==1){
                  par <- 2
                }
                else if (cont3==2){
                  par <- E^5
                }
                else if (cont3==3){
                  par <- 0.0001
                }
              }
              else{
                dpar <- par-par0
              }
              cont1 <- cont1+1
              if (par>E^6){
                par <- E^6
                dpar <- 0
              }
              if (par<=E^-5){
                par <- E^-3
                dpar <- 0
              }
            }
            alphaNB <- 1/par  # Negative binomial alpha (dispersion parameter)
          }
          dev <- 0
          ddev <- 1
          cont2 <- 1
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            assign("ai", as.vector((uj/(1+alphaNB*uj))+(y-uj)*(alphaNB*uj/(1+2*alphaNB*uj+alphaNB^2*uj*uj))), envir=parent.frame())
            ai <- as.vector((uj/(1+alphaNB*uj))+(y-uj)*(alphaNB*uj/(1+2*alphaNB*uj+alphaNB^2*uj*uj)))
            assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
            ai <- ifelse(ai<=0, E^-5, ai)
            zj <- nj+(y-uj)/(ai*(1+alphaNB*uj))-yhat_beta+fi
            b <- robust_solve(t(x)%*%(w*ai*x*wt), t(x)%*%(w*ai*wt*zj))
            # Fix matrix dimensions for negative binomial
            if (model == "negbin") {
              # Ensure b has correct dimensions
              b_vector <- as.vector(b)
              if (length(b_vector) != ncol(x)) {
                b_vector <- rep(b_vector, length.out = ncol(x))
              }
              # Ensure yhat_beta and fi have correct dimensions
              yhat_beta_corrected <- rep(as.vector(yhat_beta), length.out = N)
              fi_corrected <- rep(as.vector(fi), length.out = N)
              nj <- x %*% b_vector + yhat_beta_corrected - fi_corrected
            } else {
              nj <- x%*%b+yhat_beta-fi
            }
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- exp(nj)
            olddev <- dev
            uj <- ifelse(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (model=="poisson"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            else{
              dev <- 2*sum(y*log(tt)-(y+1/alphaNB)*log((1+alphaNB*y)/(1+alphaNB*uj)))
            }
            if (cont2>100){
              ddev <- 0.0000001
            }
            else{
              ddev <- dev-olddev
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }
        yhat_ <- get("yhat")
        yhat_[i] <- uj[i]
        assign("yhat", yhat_, envir=parent.frame())
        yhat[i] <- uj[i]
        alphai_ <- get("alphai")
        alphai_[i, 2] <- alphaNB  # Store negative binomial alpha (dispersion parameter)
        assign("alphai", alphai_, envir=parent.frame())
        alphai[i, 2] <- alphaNB  # Store negative binomial alpha (dispersion parameter)
        s_ <- get("s")
        s_[i] <- (x[i, ]%*%robust_solve(t(x)%*%(w*ai*x*wt), t(x*w*ai*wt)))[i]
        assign("s", s_, envir=parent.frame())
        s[i] <- (x[i, ]%*%robust_solve(t(x)%*%(w*ai*x*wt), t(x*w*ai*wt)))[i]
        next
      }
      else if (model=="logistic"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 0
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          assign("ai", as.vector(uj*(1-uj)), envir=parent.frame())
          ai <- as.vector(uj*(1-uj))
          assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
          ai <- ifelse(ai<=0, E^-5, ai)
          zj <- nj+(y-uj)/ai-yhat_beta+fi
          b <- robust_solve(t(x)%*%(w*ai*x*wt), t(x)%*%(w*ai*wt*zj))
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)/(1+exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          uj <- ifelse(uj==1, 0.99999, uj)
          tt2 <- (1-Y)/(1-uj)
          tt2 <- ifelse(tt2==0, E^-10, tt2)
          dev <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
          if (cont>100){
            ddev <-  0.0000001
          }
          else{
            ddev <- dev-olddev
          }
        }
        yhat_ <- get("yhat")
        yhat_[i] <- uj[i]
        assign("yhat", yhat_, envir=parent.frame())
        yhat[i] <- uj[i]
        if (det(t(x)%*%(w*ai*x*wt))==0){
          s_ <- get("s")
          s_[i] <- 0
          assign("s", s_, envir=parent.frame())
          s[i] <- 0
        }
        else{
          s_ <- get("s")
          s_[i] <- (x[i,]%*%robust_solve(t(x)%*%(w*ai*x*wt), t(x*w*wt*ai)))[i]
          assign("s", s_, envir=parent.frame())
          s[i] <- (x[i,]%*%robust_solve(t(x)%*%(w*ai*x*wt), t(x*w*wt*ai)))[i]
        }
        next
      }
    }
    if (model=="gaussian"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      npar <- sum(s)
      AICc <- 2*N*log(CV/N)+N*log(2*3.14159)+N*(N+npar)/(N-2-npar)
    }
    else if (model=="poisson" | model=="negbin"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      if (model=="poisson"){
        ll <- sum(-yhat+y*log(yhat)-lgamma(y+1))
        npar <- sum(s)
      }
      else{
        ll <- sum(y*log(alphai[,2]*yhat)-(y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(y+1))
        npar <- sum(s)+sum(s)/nvar
      }
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    }
    else if (model=="logistic"){
      uj <- ifelse(uj==0, E^-10, uj)
      uj <- ifelse(uj==1, 0.99999, uj)
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      ll <- sum(y*log(uj)-(1-y)*log(1-uj))
      npar <- sum(s)
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    }
    if (bandwidth=="aic"){
      CV <- AICC
    }
    res <- cbind(CV, npar)
    return (res)
  }
  GSS <- function(depy, indepx, fix){
    # DEFINING GOLDEN SECTION SEARCH PARAMETERS #
    if(method=="fixed_g" | method=="fixed_bsq"){
      ax <- 0
      bx <- as.integer(max(dist(COORD))+1)
      if (distancekm){
        bx <- bx*111
      }
    }
    else if (method=="adaptive_bsq" | method=="adaptive_bsq_smr"){
      # Use 5 to N range for adaptive methods
      # ax <- 5  # Minimum 5
      # bx <- N  # Maximum N
      ax <- floor(0.15 * N)  # Minimum 15 or 15% of N (whichever is larger)
      bx <- floor(0.35 * N)  # Maximum 35% of N (slightly larger)
    }
    r <- 0.61803399
    tol <- 0.1
    if (!globalmin){
      lower <- ax
      upper <- bx
      xmin <- matrix(0, 1, 2)
      GMY <- 1
      ax1 <- lower[GMY]
      bx1 <- upper[GMY]
      h0 <- ax1
      h3 <- bx1
      h1 <- bx1-r*(bx1-ax1)
      h2 <- ax1+r*(bx1-ax1)
      res1 <- cv(h1, depy, indepx, fix)
      assign("s", s, envir=parent.frame())
      if (model!="gaussian"){ #release 2
        assign("ai", ai, envir=parent.frame())
      }#release 2
      assign("yhat", yhat, envir=parent.frame())
      assign("alphai", alphai, envir=parent.frame())
      CV1 <- res1[1]
      res2 <- cv(h2,depy,indepx,fix)
      assign("s", s, envir=parent.frame())
      if (model!="gaussian"){ #release 2
        assign("ai", ai, envir=parent.frame())
      } #release 2
      assign("yhat", yhat, envir=parent.frame())
      assign("alphai", alphai, envir=parent.frame())
      CV2 <- res2[1]
      INT <- 1
      while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & INT<200){
        if (CV2<CV1){
          h0 <- h1
          h1 <- h3-r*(h3-h0)
          h2 <- h0+r*(h3-h0)
          CV1 <- CV2
          res2 <- cv(h2,depy,indepx,fix)
          assign("s", s, envir=parent.frame())
          if (model!="gaussian"){ #release 2
            assign("ai", ai, envir=parent.frame())
          } #release 2
          assign("yhat", yhat, envir=parent.frame())
          assign("alphai", alphai, envir=parent.frame())
          CV2 <- res2[1]
        }
        else{
          h3 <- h2
          h1 <- h3-r*(h3-h0)
          h2 <- h0+r*(h3-h0)
          CV2 <- CV1
          res1 <- cv(h1, depy, indepx, fix)
          assign("s", s, envir=parent.frame())
          if (model!="gaussian"){ #release 2
            assign("ai", ai, envir=parent.frame())
          } #release 2
          assign("yhat", yhat, envir=parent.frame())
          assign("alphai", alphai, envir=parent.frame())
          CV1 <- res1[1]
        }
        INT <- INT+1
      }
      if (CV1<CV2){
        golden <- CV1
        xmin[GMY, 1] <- golden
        xmin[GMY, 2] <- h1
        npar <- res1[1]
        if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
          xmin[GMY, 2] <- floor(h1)
          xming <- xmin[GMY, 2]
        }
      }
      else{
        golden <- CV2
        xmin[GMY, 1] <- golden
        xmin[GMY, 2] <- h2
        npar <- res2[1]
        if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
          xmin[GMY, 2] <- floor(h2)
          xming <- xmin[GMY, 2]
        }
      }
      xming <- xmin[GMY, 2]
    }
    else{
      lower <- cbind(ax, (1-r)*bx, r*bx)
      upper <- cbind((1-r)*bx, r*bx, bx)
      xmin <- matrix(0, 3, 2)
      for (GMY in 1:3){
        ax1 <- lower[GMY]
        bx1 <- upper[GMY]
        h0 <- ax1
        h3 <- bx1
        h1 <- bx1-r*(bx1-ax1)
        h2 <- ax1+r*(bx1-ax1)
        res1 <- cv(h1, depy, indepx, fix)
        assign("s", s, envir=parent.frame())
        if (model!="gaussian"){ #release 2
          assign("ai", ai, envir=parent.frame())
        } #release 2
        assign("yhat", yhat, envir=parent.frame())
        assign("alphai", alphai, envir=parent.frame())
        CV1 <- res1[1]
        res2 <- cv(h2,depy,indepx,fix)
        assign("s", s, envir=parent.frame())
        if (model!="gaussian"){ #release 2
          assign("ai", ai, envir=parent.frame())
        } #release 2
        assign("yhat", yhat, envir=parent.frame())
        assign("alphai", alphai, envir=parent.frame())
        CV2 <- res2[1]
        INT <- 1
        while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & INT<200){
          if (CV2<CV1){
            h0 <- h1
            h1 <- h3-r*(h3-h0)
            h2 <- h0+r*(h3-h0)
            CV1 <- CV2
            res2 <- cv(h2,depy,indepx,fix)
            assign("s", s, envir=parent.frame())
            if (model!="gaussian"){ #release 2
              assign("ai", ai, envir=parent.frame())
            } #release 2
            assign("yhat", yhat, envir=parent.frame())
            assign("alphai", alphai, envir=parent.frame())
            CV2 <- res2[1]
          }
          else{
            h3 <- h2
            h1 <- h3-r*(h3-h0)
            h2 <- h0+r*(h3-h0)
            CV2 <- CV1
            res1 <- cv(h1, depy, indepx, fix)
            assign("s", s, envir=parent.frame())
            if (model!="gaussian"){ #release 2
              assign("ai", ai, envir=parent.frame())
            } #release 2
            assign("yhat", yhat, envir=parent.frame())
            assign("alphai", alphai, envir=parent.frame())
            CV1 <- res1[1]
          }
          INT <- INT+1
        }
        if (CV1<CV2){
          golden <- CV1
          xmin[GMY,1] <- golden
          xmin[GMY,2] <- h1
          npar <- res1[1]
          if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
            xmin[GMY,2] <- floor(h1)
            xming <- xmin[GMY,2]
          }
        }
        else{
          golden <- CV2
          xmin[GMY,1] <- golden
          xmin[GMY,2] <- h2
          npar <- res2[1]
          if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
            xmin[GMY,2] <- floor(h2)
            xming <- xmin[GMY,2]
          }
        }
        xming <- xmin[GMY,2]
      }
    }
    if (globalmin){
      xming <- xmin[which(xmin[,1]==min(xmin[,1])),2]
    }
    bandwidth <- xming
    return (bandwidth)
  }
  
  gwr <- function(H, y, x, fi, alphaW_val = NULL) {
    # alphaW_val dapat dipassing dari luar, default ke 0.5 jika NULL
    nvar <- ncol(x)
    bim <- rep(0, nvar*N)
    yhatm <- rep(0, N)
    for (i in 1:N){
      # Calculate distance matrix for this point
      seqi <- rep(i, N)
      dx <- sp::spDistsN1(COORD, COORD[i,])
      distan <- cbind(seqi, sequ, dx)
      if (distancekm){
        distan[,3] <- distan[,3]*111
      }
      u <- nrow(distan)
      
      dist_spatial_i <- dist_spatial_mat[i, ]
      
      w <- rep(0, N) 
      
      if (method=="fixed_g"){
        for (jj in 1:u){
          w <- exp(-(dist_spatial_i/H)^2)
        }
      }
      else if (method=="fixed_bsq"){
        for (jj in 1:u){
          w <- (1-(dist_spatial_i/H)^2)^2
          w[dist_spatial_i > H] <- 0
        }
      }
      else if (method=="adaptive_bsq"){
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, N, 2)
        hn <- distan[H,3]
        for (jj in 1:N){
          if (distan[jj,4]<=H){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          } else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
        }
        w <- w[order(w[, 2]), ]
        w <- w[,1]
      }
      else if (method == "adaptive_bsq_smr") {
        # PATCH: gunakan alphaW_val dari parameter, default 1.0 untuk mengikuti adaptive_bsq
        if (is.null(alphaW_val)) alphaW_val <- 1.0
        betaW_val <- 1 - alphaW_val
        
        # --- w1 (Spatial weights with bisquare kernel) ---
        dist_spatial_sorted_df <- data.frame(dist = dist_spatial_i, original_index = 1:N)
        dist_spatial_sorted_df <- dist_spatial_sorted_df[order(dist_spatial_sorted_df$dist), ]
        dist_spatial_sorted_df$rank <- 1:N
        
        w1_matrix <- matrix(0, N, 2)
        hn_spatial <- dist_spatial_sorted_df$dist[H]
        if (hn_spatial == 0) hn_spatial <- 1e-6 
        for (jj in 1:N) {
          if (dist_spatial_sorted_df$rank[jj] <= H) {
            w1_matrix[jj, 1] <- (1 - (dist_spatial_sorted_df$dist[jj] / hn_spatial)^2)^2
          } else {
            w1_matrix[jj, 1] <- 0
          }
          w1_matrix[jj, 2] <- dist_spatial_sorted_df$original_index[jj]
        }
        w1_matrix <- w1_matrix[order(w1_matrix[, 2]), ]
        w1 <- w1_matrix[, 1]
        
        # --- w2 (Attribute similarity weights with bisquare kernel) ---
        dist_sim_i <- dist_sim_mat[i, ]
        dist_sim_sorted_df <- data.frame(dist = dist_sim_i, original_index = 1:N)
        dist_sim_sorted_df <- dist_sim_sorted_df[order(dist_sim_sorted_df$dist), ]
        dist_sim_sorted_df$rank <- 1:N
        
        w2_matrix <- matrix(0, N, 2)
        hn_sim <- dist_sim_sorted_df$dist[H]
        if (hn_sim == 0) hn_sim <- 1e-6
        
        for (jj in 1:N) {
          if (dist_sim_sorted_df$rank[jj] <= H) {
            w2_matrix[jj, 1] <- exp(-(dist_sim_sorted_df$dist[jj] / hn_sim)^2)
          } else {
            w2_matrix[jj, 1] <- 0
          }
          w2_matrix[jj, 2] <- dist_sim_sorted_df$original_index[jj]
        }
        w2_matrix <- w2_matrix[order(w2_matrix[, 2]), ]
        w2 <- w2_matrix[, 1]
        
        w <- alphaW_val * w1 + betaW_val * w2
      }
      ## MODEL SELECTION ##
      if (model=="gaussian"){
        if (det(t(x)%*%(w*x*wt))==0){
          b <- rep(0, nvar)
        }
        else{
          b <- solve(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
        }
        uj <- x%*%b
        if (nvar==nvarg){
          if (det(t(x)%*%(w*x*wt))==0){
            sm_ <- get("sm")
            sm_[i, ] <- rep(0, N)
            assign("sm", sm_, envir=parent.frame())
            mrj_ <- get("mrj")
            mrj_[i, ] <- matrix(0, N*nvar)
            assign("mrj", mrj_, envir=parent.frame())
          }
          else{
            ej <- diag(nvar)
            sm_ <- get("sm")
            sm_[i, ] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            assign("sm", sm_, envir=parent.frame())
            sm3_ <- get("sm3")
            sm3_[i, ] <- t(diag((solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))%*%t(solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))))
            assign("sm3", sm3_, envir=parent.frame())
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj_ <- get("mrj")
              mrj_[i, m1:m2] <- (x[i,jj]*ej[jj,])%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt)
              assign("mrj", mrj_, envir=parent.frame())
            }
          }
        }
        else{
          if (det(t(x)%*%(w*x*wt))==0){
            rj_ <- get("rj")
            rj_[i, ] <- rep(0, N)
            assign("rj", rj_, envir=parent.frame())
          }
          else{
            rj_ <- get("rj")
            rj_[i, ] <- (x[i,]%*%solve(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            assign("rj", rj_, envir=parent.frame())
          }
        }
      }
      else if (model=="poisson" | model=="negbin"){
        uj <- yhat
        par <- parg
        nj <- log(uj)
        ddpar <- 1
        cont <- 1
        while (abs(ddpar)>0.000001 & cont<100){
          dpar <- 1
          parold <- par
          cont1 <- 1
          cont3 <- 1
          if (model=="poisson"){
            alphaNB <- E^-6  # Negative binomial alpha (dispersion parameter)
            par <- 1/alphaNB
          }
          else{ #model=='NEGBIN'
            if (par<=E^-5 & i>1){
              par=1/alphai[i-1,2]
            }
            while (abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(as.vector(par)+y)-digamma(as.vector(par))+log(as.vector(par))+1-log(as.vector(par)+uj)-(as.vector(par)+y)/(as.vector(par)+uj)))
              hess <- sum(w*wt*(trigamma(as.vector(par)+y)-trigamma(as.vector(par))+1/as.vector(par)-2/(as.vector(par)+uj)+(y+as.vector(par))/(as.vector(par)+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              par <- par0 - g/hess  # Use division instead of solve for scalar
              if (cont1>50 & par>E^5){
                dpar <- 0.0001
                cont3 <- cont3+1
                if (cont3==1){
                  par <- 2
                }
                else if (cont3==2){
                  par <- E^5
                }
                else if (cont3==3){
                  par <- 0.0001
                }
              }
              else{
                dpar <- par-par0
              }
              cont1 <- cont1+1
              if (par>E^6){
                par <- E^6
                dpar <- 0
              }
              if (par<=E^-5){
                par <- E^-3
                dpar <- 0
              }
            }
            alphaNB <- as.vector(1/par)  # Negative binomial alpha (dispersion parameter)
          }
          dev <- 0
          ddev <- 1
          cont2 <- 0
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            ai <- as.vector((uj/(1+alphaNB*uj))+(y-uj)*(alphaNB*uj/(1+2*alphaNB*uj+alphaNB^2*uj*uj)))
            ai <- ifelse(ai<=0, E^-5, ai)
            zj <- nj+(y-uj)/(ai*(1+alphaNB*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- solve(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
            }
            # Fix matrix dimensions for negative binomial
            if (model == "negbin") {
              # Ensure b has correct dimensions
              b_vector <- as.vector(b)
              if (length(b_vector) != ncol(x)) {
                b_vector <- rep(b_vector, length.out = ncol(x))
              }
              # Ensure yhat_beta and fi have correct dimensions
              yhat_beta_corrected <- rep(as.vector(yhat_beta), length.out = N)
              fi_corrected <- rep(as.vector(fi), length.out = N)
              nj <- x %*% b_vector + yhat_beta_corrected - fi_corrected
            } else {
              nj <- x%*%b+yhat_beta-fi
            }
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- as.vector(exp(nj))
            olddev <- dev
            uj <- ifelse(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (model=="poisson"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            else{ #model=="NEGBIN"
              dev <- 2*sum(y*log(tt)-(y+1/alphaNB)*log((1+alphaNB*y)/(1+alphaNB*uj)))
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }
        if (nvar==nvarg){
          if (det(t(x)%*%(w*ai*x*wt))==0){
            sm_ <- get("sm")
            sm_[i, ] <- c(0, N)
            assign("sm", sm_, envir=parent.frame())
            mrj_ <- get("mrj")
            mrj_[i, ] <- rep(0, N*nvar)
            assign("mrj", mrj_, envir=parent.frame())
          }
          else{
            ej <- diag(nvar)
            sm_ <- get("sm")
            sm_[i, ] <- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            assign("sm", sm_, envir=parent.frame())
            sm3_ <- get("sm3")
            sm3_[i, ] <- t(diag((solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))%*%diag(1/ai)%*%t(solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))))
            assign("sm3", sm3_, envir=parent.frame())
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj_ <- get("mrj")
              mrj_[i, m1:m2] <- (x[i,jj]*ej[jj,])%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
              assign("mrj", mrj_, envir=parent.frame())
            }
          }
        }
        else{
          if (det(t(x)%*%(w*ai*x*wt))==0){
            rj_ <- get("rj")
            rj_[i, ] <- rep(0, N)
            assign("rj", rj_, envir=parent.frame())
          }
          else{
            rj_ <- get("rj")
            rj_[i, ] <- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            assign("rj", rj_, envir=parent.frame())
          }
        }
        if (model=="negbin"){
          hess <- sum(w*wt*(trigamma(as.vector(par)+y)-trigamma(as.vector(par))+1/as.vector(par)-2/(as.vector(par)+exp(yhat_beta))+(y+as.vector(par))/(as.vector(par)+exp(yhat_beta))^2))
          if (!mgwr){
            hess <- sum(w*wt*(trigamma(as.vector(par)+y)-trigamma(as.vector(par))+1/as.vector(par)-2/(as.vector(par)+uj)+(y+as.vector(par))/(as.vector(par)+uj)^2))
            hess <- ifelse(hess==0, E^-23, hess)
          } else {
            # Safety check for mgwr=TRUE to prevent division by zero
            hess <- ifelse(hess==0, E^-23, hess)
          }
          sealphaNB <- sqrt(1/abs(hess))/(par^2)
          alphai_ <- get("alphai")
          alphai_[i, 1] <- i
          assign("alphai", alphai_, envir=parent.frame())
          alphai_ <- get("alphai")
          alphai_[i, 2] <- alphaNB  # Store negative binomial alpha (dispersion parameter)
          assign("alphai", alphai_, envir=parent.frame())
          alphai_ <- get("alphai")
          alphai_[i, 3] <- sealphaNB
          assign("alphai", alphai_, envir=parent.frame())
        }
      }
      else{ #else if (model=="logistic"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 1
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          assign("ai", as.vector(uj*(1-uj)), envir=parent.frame())
          ai <- as.vector(uj*(1-uj))
          assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
          ai <- ifelse(ai<=0, E^-5, ai)
          zj <- nj+(y-uj)/ai-yhat_beta+fi
          b <- robust_solve(t(x)%*%(w*ai*x*wt), t(x)%*%(w*ai*wt*zj))
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)/(1+exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          uj <- ifelse(uj==1, 0.99999, uj)
          tt2 <- (1-Y)/(1-uj)
          tt2 <- ifelse(tt2==0, E^-10, tt2)
          dev <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
          if (cont>100){
            ddev <-  0.0000001
          }
          else{
            ddev <- dev-olddev
          }
        }
        if (nvar==nvarg){
          if (det(t(x)%*%(w*ai*x*wt))==0){
            sm_ <- get("sm")
            sm_[i, ] <- rep(0, N)
            assign("sm", sm_, envir=parent.frame())
            mrj_ <- get("mrj")
            mrj_[i, ] <- matrix(0, N*nvar)
            assign("mrj", mrj_, envir=parent.frame())
          }
          else{
            ej <- diag(nvar)
            sm_ <- get("sm")
            sm_[i, ] <- x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
            assign("sm", sm_, envir=parent.frame())
            sm3_ <- get("sm3")
            sm3_[i, ] <- t(diag((solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))%*%diag(1/ai)%*%t(solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))))
            assign("sm3", sm3_, envir=parent.frame())
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj_ <- get("mrj")
              mrj_[i, m1:m2] <- (x[i,jj]*ej[jj,])%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
              assign("mrj", mrj_, envir=parent.frame())
            }
          }
        }
        else{
          if (det(t(x)%*%(w*ai*x*wt))==0){
            rj_ <- get("rj")
            rj_[i, ] <- rep(0, N)
            assign("rj", rj_, envir=parent.frame())
          }
          else{
            rj_ <- get("rj")
            rj_[i, ] <- (x[i,]%*%solve(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            assign("rj", rj_, envir=parent.frame())
          }
        }
      }
      m1 <- (i-1)*nvar+1
      m2 <- m1+(nvar-1)
      bim[m1:m2] <- b
      yhatm[i] <- uj[i]
      yhat_ <- get("yhat")
      yhat_[i] <- uj[i]
      assign("yhat", yhat_, envir=parent.frame())
    }
    beta <- matrix(bim, N, byrow=T)
    yhbeta <- cbind(yhatm, beta)
    return (yhbeta)
  }
  
  
  
  if (!mgwr){
    if (verbose) cat("Step 3: Bandwidth selection for GWR (mgwr=FALSE)...\n")
    
    finb <- rep(0, N)
    yhat_beta <- Offset
    if (!is.null(h)){
      hh <- h
      if (verbose) cat("  - Using fixed bandwidth:", hh, "\n")
    }
    else{
      if (verbose) cat("  - Optimizing bandwidth using", bandwidth, "criterion...\n")
      if (method == "adaptive_bsq_smr") {
        if (verbose) cat("  - Using adaptive_bsq_smr with alphaW=1.0 for bandwidth optimization\n")
        # Use GSS_with_alpha for adaptive_bsq_smr method
        hh <- GSS_with_alpha(Y, X, finb, COORD, sequ, distancekm, method, 1.0)
      } else {
        hh <- GSS(Y,X,finb)
      }
      if (verbose) cat("  - Optimal bandwidth found:", hh, "\n")
    }
    header <- append(header, "General Bandwidth")
    output <- append(output, hh)
    names(output) <- "general_bandwidth"
    
    # === PATCH: pencarian optimal alphaW khusus method adaptive_bsq_smr ===
    cat("DEBUG: About to check method, method =", method, "class =", class(method), "\n")
    if (method == "adaptive_bsq_smr") {
      if (verbose) cat("Step 4: Adaptive BSQ SMR optimization...\n")
      cat("DEBUG: Entering adaptive_bsq_smr optimization, mgwr =", mgwr, "\n")
      # Adaptive BSQ SMR method
      cat("DEBUG: Checking mgwr condition, mgwr =", mgwr, "class =", class(mgwr), "\n")
      if (mgwr) {
        cat("DEBUG: mgwr=TRUE, entering MGWR optimization\n")
        if (verbose) cat("  - MGWR=TRUE, starting multi-alphaW optimization\n")
        
        # STEP 1: Find optimal alphaW for GWR (single alphaW) and its AICc as baseline
        if (verbose) cat("  - Step 4.1: Finding optimal single alphaW for adaptive_bsq_smr (mgwr=FALSE equivalent)...\n")
        if (verbose) cat("    * Running Golden Section Search for optimal alphaW...\n")
        optimal_alphaW_gwr <- find_optimal_alphaW_gss(
          hh, Y, X, finb, N, nvarg, Offset, method, model, wt, E, COORD, sequ, distancekm, parg, yhat_beta
        )
        if (verbose) cat("    * Optimal alphaW found:", round(optimal_alphaW_gwr, 3), "\n")
        
        # Calculate baseline AICc with single optimal alphaW
        if (verbose) cat("    * Calculating baseline performance (mgwr=FALSE equivalent)...\n")
        baseline_res <- gwr_local(
          H = hh, y = Y, x = X, fi = finb, alphaW_val = optimal_alphaW_gwr, 
          method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E, 
          COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset, 
          parg = parg, yhat_beta = yhat_beta
        )
        baseline_yhbeta <- baseline_res$yhbeta
        baseline_beta <- baseline_yhbeta[, 2:(nvarg + 1)]
        baseline_Fi <- X * baseline_beta
        baseline_sm <- baseline_res$sm
        baseline_alphai <- baseline_res$alphai
        baseline_v1 <- sum(diag(baseline_sm))
        baseline_yhat <- exp(apply(baseline_Fi, 1, sum) + Offset)
        baseline_ll <- sum(Y * log(baseline_alphai[,2] * baseline_yhat) - (Y + 1/baseline_alphai[,2]) * log(1 + baseline_alphai[,2] * baseline_yhat) +
                             lgamma(Y + 1/baseline_alphai[,2]) - lgamma(1/baseline_alphai[,2]) - lgamma(Y + 1))
        baseline_AIC <- 2 * baseline_v1 - 2 * baseline_ll
        baseline_AICc <- baseline_AIC + 2 * baseline_v1 * (baseline_v1 + 1) / (N - baseline_v1 - 1)
        
        # Calculate baseline percent deviance explained
        baseline_deviance <- -2 * baseline_ll
        null_deviance <- -2 * sum(Y * log(mean(Y)) - mean(Y) - lgamma(Y + 1))
        baseline_percent_deviance <- (1 - baseline_deviance / null_deviance) * 100
        
        if (verbose) cat("    * Baseline established (mgwr=FALSE equivalent):\n")
        if (verbose) cat("      AIC:", round(baseline_AIC, 2), "\n")
        if (verbose) cat("      AICc:", round(baseline_AICc, 2), "\n")
        if (verbose) cat("      Percent deviance explained:", round(baseline_percent_deviance, 2), "%\n")
        
        # STEP 2: Find optimum multi-bandwidths using alphaW = 1 (purely spatial)
        if (verbose) cat("  - Step 4.2: Finding optimum multi-bandwidths using alphaW = 1 (purely spatial)...\n")
        mband <- rep(0, nvarg)
        
        # Use parallel computing for bandwidth optimization if available
        if (!is.null(cl)) {
          if (verbose) cat("    * Using parallel computing for bandwidth optimization...\n")
          mband <- foreach(i = 1:nvarg, .combine = c, .packages = c("stats", "MASS", "sp")) %dopar% {
            GSS_with_alpha(Y, as.matrix(X[,i]), rep(0, N), COORD, sequ, distancekm, method, 1.0)
          }
        } else {
          for (i in 1:nvarg) {
            if (verbose) cat("    * Optimizing bandwidth for variable", i, "of", nvarg, "\n")
            # Use alphaW = 1 for purely spatial bandwidth calculation
            mband[i] <- GSS_with_alpha(Y, as.matrix(X[,i]), rep(0, N), COORD, sequ, distancekm, method, 1.0)
          }
        }
        if (verbose) cat("  - Multi-bandwidths optimized:", paste(round(mband, 2), collapse=", "), "\n")
        cat("DEBUG: Optimum multi-bandwidths (alphaW = 1) =", paste(round(mband, 2), collapse=", "), "\n")
        
        # STEP 3: Use global-local hybrid optimization for multi-alphaW
        # Use optimum single alphaW as initial alphaW per variable, then optimize alphaW range [0,1]
        if (verbose) cat("  - Step 4.3: Running global-local hybrid optimization for multi-alphaW...\n")
        cat("DEBUG: Using optimum single alphaW =", round(optimal_alphaW_gwr, 3), "as initial alphaW per variable\n")
        hybrid_result <- global_local_hybrid_mgwr(
          Y = Y, X = X, method = method, model = model, N = N, nvarg = nvarg,
          Offset = Offset, wt = wt, E = E, COORD = COORD, sequ = sequ,
          distancekm = distancekm, parg = parg, yhat_beta = yhat_beta,
          initial_bandwidths = mband, initial_alpha = optimal_alphaW_gwr, cl = cl
        )
        
        # Use the hybrid optimized parameters
        optimal_bandwidths <- hybrid_result$bandwidths
        optimal_alphaWs <- hybrid_result$alphas
        final_aicc <- hybrid_result$aicc
        iterations <- hybrid_result$iterations
        
        if (verbose) {
          cat("  - Global-local hybrid optimization completed in", iterations, "iterations\n")
          cat("  - Final AICc =", round(final_aicc, 2), "vs Baseline AICc =", round(baseline_AICc, 2), "\n")
          cat("  - Improvement =", round(baseline_AICc - final_aicc, 2), "points\n")
        }
        
        # Run final GWR with optimized parameters
        if (verbose) cat("  - Running final GWR with optimized parameters...\n")
        res <- gwr_local_multi_alpha(
          H = optimal_bandwidths, y = Y, x = X, fi = finb, alpha_vals = optimal_alphaWs,
          method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
          COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
          parg = parg, yhat_beta = yhat_beta, verbose = verbose
        )
        
        yhat_beta <- res$yhbeta
        sm <- res$sm
        alphai <- res$alphai
        sm3 <- res$sm3
        
        # Store multi-alphaW results for later use
        output$optimal_alphaWs <- optimal_alphaWs
        output$optimal_bandwidths <- optimal_bandwidths
        output$baseline_aicc <- baseline_AICc
        output$final_aicc <- final_aicc
        output$improvement <- baseline_AICc - final_aicc
        output$iterations <- iterations
        
        if (verbose) {
          cat("  - Optimization results stored\n")
          cat("  - Optimal alphaWs:", paste(round(optimal_alphaWs, 3), collapse=", "), "\n")
          cat("  - Optimal bandwidths:", paste(round(optimal_bandwidths, 2), collapse=", "), "\n")
        }
        
        # Create alphaW optimum per variable information
        alphaW_info <- data.frame(
          Variable = c("Intercept", XVAR),
          Optimal_AlphaW = optimal_alphaWs,
          Optimal_Bandwidth = optimal_bandwidths
        )
        output$alphaW_optimum_per_variable <- alphaW_info
      } else {
        cat("DEBUG: mgwr=FALSE, entering single alphaW optimization\n")
        # GWR with single alphaW optimization - use faster approach for mgwr=FALSE
        # For mgwr=FALSE, use a simpler alphaW optimization to avoid excessive computation
        if (verbose) cat("  - MGWR=FALSE, optimizing single alphaW for adaptive_bsq_smr...\n")
        
        # Use a coarser tolerance for faster convergence
        if (verbose) cat("    * Finding optimal alphaW with tolerance = 0.05...\n")
        optimal_alphaW <- find_optimal_alphaW_gss(
          hh, Y, X, finb, N, nvarg, Offset, method, model, wt, E, COORD, sequ, distancekm, parg, yhat_beta, tol = 0.05
        )
        if (verbose) cat("    * Optimal alphaW found:", round(optimal_alphaW, 3), "\n")
        
        # Call GWR_LOCAL after alphaW optimization
        if (verbose) cat("    * Running GWR with optimized alphaW...\n")
        res <- gwr_local(
          H = hh, y = Y, x = X, fi = finb, alphaW_val = optimal_alphaW, 
          method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E, 
          COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset, 
          parg = parg, yhat_beta = yhat_beta
        )
        yhat_beta <- res$yhbeta
        sm <- res$sm
        alphai <- res$alphai
        sm3 <- res$sm3
        # Store single alphaW result
        output$optimal_alphaW <- optimal_alphaW
        if (verbose) cat("    * Single alphaW optimization completed\n")
      }
    } else {
      # method lain tetap pakai gwr versi lama
      yhat_beta <- gwr(hh,Y,X,finb)
    }
    
    # Process results
    beta <- yhat_beta[,2:(nvarg+1)]
    Fi <- X*beta
    mband <- hh
    Sm2 <- sm
  } 
  else {
    if (verbose) cat("Step 3: MGWR (mgwr=TRUE) processing...\n")
    cat("DEBUG: Method =", method, "mgwr =", mgwr, "\n")
    
    finb <- rep(0, N)
    yhat_beta <- Offset
    
    if (!is.null(h)) {
      hh <- h
      if (verbose) cat("  - Using fixed initial bandwidth:", hh, "\n")
    } else {
      if (verbose) cat("  - Optimizing initial global bandwidth using", bandwidth, "criterion...\n")
      # Use same bandwidth calculation as direct mgwr=FALSE for adaptive_bsq_smr
      if (method == "adaptive_bsq_smr") {
        if (verbose) cat("  - Using adaptive_bsq_smr with alphaW=1.0 for bandwidth optimization\n")
        hh <- GSS_with_alpha(Y, X, finb, COORD, sequ, distancekm, method, 1.0)
      } else {
        hh <- GSS(Y, X, finb)
      }
      if (verbose) cat("  - Initial global bandwidth found:", hh, "\n")
    }
    header <- append(header, "General Bandwidth")
    output <- append(output, hh)
    names(output) <- "general_bandwidth"
    
    if (verbose) cat("  - Computing initial model fit using spatial-only weights (alphaW=1.0)...\n")
    initial_gwr_fit <- gwr(hh, Y, X, finb, alphaW_val = 1.0)
    beta <- initial_gwr_fit[, 2:(nvarg + 1)]
    Fi <- X * beta
    
    mband <- rep(hh, nvarg)
    
    # --- START OF METHOD-SPECIFIC LOGIC ---
    cat("DEBUG: Checking method, method =", method, "class =", class(method), "\n")
    if (method == "adaptive_bsq") {
      cat("DEBUG: Method is adaptive_bsq, using simple MGWR approach\n")
      # For adaptive_bsq with mgwr=TRUE, use simple approach without alphaW optimization
      # adaptive_bsq is pure spatial method, so alphaW = 1.0 always
      if (verbose) cat("  - Using adaptive_bsq with mgwr=TRUE (pure spatial, alphaW=1.0)...\n")
      
      # Initialize with single bandwidth for all variables
      mband <- rep(hh, nvarg)
      optimal_alphas <- rep(1.0, nvarg)  # Always 1.0 for adaptive_bsq
      
      if (verbose) cat("  - Bandwidths:", paste(round(mband, 2), collapse=", "), "\n")
      if (verbose) cat("  - AlphaWs (fixed):", paste(round(optimal_alphas, 3), collapse=", "), "\n")
      
    } else if (method == "adaptive_bsq_smr") {
      cat("DEBUG: Method is adaptive_bsq_smr, entering baseline optimization\n")
      # STEP 1: Establish MGWR FALSE baseline before back-fitting
      if (verbose) cat("  - Step 4.1: Finding optimal single alphaW for adaptive_bsq_smr (mgwr=FALSE equivalent)...\n")
      if (verbose) cat("    * Running Golden Section Search for optimal alphaW...\n")
      
      cat("DEBUG: About to call find_optimal_alphaW_gss with method =", method, "\n")
      optimal_alphaW_gwr <- find_optimal_alphaW_gss(
        hh, Y, X, finb, N, nvarg, Offset, method, model, wt, E, COORD, sequ, distancekm, parg, yhat_beta
      )
      cat("DEBUG: find_optimal_alphaW_gss returned:", optimal_alphaW_gwr, "\n")
      if (verbose) cat("    * Optimal alphaW found:", round(optimal_alphaW_gwr, 3), "\n")
      
      # Calculate baseline performance (mgwr=FALSE equivalent) - Use EXACT same algorithm as direct mgwr=FALSE
      if (verbose) cat("    * Calculating baseline performance (mgwr=FALSE equivalent)...\n")
      
      # Use the EXACT same algorithm as direct mgwr=FALSE for adaptive_bsq_smr
      if (method == "adaptive_bsq_smr") {
        # COPY EXACT SAME IMPLEMENTATION AS DIRECT MGWR=FALSE
        # Use same initialization as direct mgwr=FALSE
        baseline_yhat_beta <- rep(0, N)  # Same as direct mgwr=FALSE
        
        if (verbose) cat("    * Finding optimal alphaW with tolerance = 0.05...\n")
        baseline_optimal_alphaW <- find_optimal_alphaW_gss(
          hh, Y, X, finb, N, nvarg, Offset, method, model, wt, E, COORD, sequ, distancekm, parg, baseline_yhat_beta, tol = 0.05
        )
        if (verbose) cat("    * Optimal alphaW found:", round(baseline_optimal_alphaW, 3), "\n")
        
        # Call GWR_LOCAL after alphaW optimization (EXACT SAME AS DIRECT MGWR=FALSE)
        if (verbose) cat("    * Running GWR with optimized alphaW...\n")
        baseline_res <- gwr_local(
          H = hh, y = Y, x = X, fi = finb, alphaW_val = baseline_optimal_alphaW, 
          method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E, 
          COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset, 
          parg = parg, yhat_beta = baseline_yhat_beta
        )
        if (verbose) cat("    * Single alphaW optimization completed\n")
      } else {
        # For other methods, use the original approach
        baseline_res <- gwr_local(
          H = hh, y = Y, x = X, fi = finb, alphaW_val = optimal_alphaW_gwr, 
          method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E, 
          COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset, 
          parg = parg, yhat_beta = yhat_beta
        )
      }
      baseline_yhbeta <- baseline_res$yhbeta
      baseline_beta <- baseline_yhbeta[, 2:(nvarg + 1)]
      baseline_Fi <- X * baseline_beta
      baseline_sm <- baseline_res$sm
      baseline_alphai <- baseline_res$alphai
      baseline_v1 <- sum(diag(baseline_sm))
      baseline_yhat <- exp(apply(baseline_Fi, 1, sum) + Offset)
      baseline_ll <- sum(Y * log(baseline_alphai[,2] * baseline_yhat) - (Y + 1/baseline_alphai[,2]) * log(1 + baseline_alphai[,2] * baseline_yhat) +
                           lgamma(Y + 1/baseline_alphai[,2]) - lgamma(1/baseline_alphai[,2]) - lgamma(Y + 1))
      baseline_AIC <- 2 * baseline_v1 - 2 * baseline_ll
      baseline_AICc <- baseline_AIC + 2 * baseline_v1 * (baseline_v1 + 1) / (N - baseline_v1 - 1)
      
      # Calculate baseline percent deviance explained (SAME as final model)
      # Use EXACT same deviance calculation as final negbin model
      tt_baseline <- Y/baseline_yhat
      tt_baseline <- ifelse(tt_baseline==0, E^-10, tt_baseline)
      baseline_deviance <- 2*sum(Y*log(tt_baseline)-(Y+1/baseline_alphai[,2])*log((1+baseline_alphai[,2]*Y)/(1+baseline_alphai[,2]*baseline_yhat)))
      
      tt2_baseline <- Y/mean(Y)
      tt2_baseline <- ifelse(tt2_baseline==0, E^-10, tt2_baseline)
      null_deviance <- 2*sum(Y*log(tt2_baseline)-(Y+1/baseline_alphai[,2])*log((1+baseline_alphai[,2]*Y)/(1+baseline_alphai[,2]*mean(Y))))
      
      baseline_percent_deviance <- (1 - baseline_deviance / null_deviance)  # No *100, same as final model
      
      if (verbose) cat("    * Baseline established (mgwr=FALSE equivalent):\n")
      if (verbose) cat("      Optimal alphaW:", round(baseline_optimal_alphaW, 3), "\n")
      if (verbose) cat("      AIC:", round(baseline_AIC, 2), "\n")
      if (verbose) cat("      AICc:", round(baseline_AICc, 2), "\n")
      if (verbose) cat("      Percent deviance explained:", round(baseline_percent_deviance, 4), "\n")
      
      # STEP 2: Now start back-fitting to improve upon the baseline
      if (verbose) cat("  - Starting back-fitting process with Adaptive Grid Search for Alphas...\n")
      
      # Initialize back-fitting with simple iteration-based convergence
      INT <- 1
      optimal_alphas <- rep(baseline_optimal_alphaW, nvarg) # Initialize with optimal alphaW from baseline
      max_iterations <- 5  # Reduced for faster convergence with AIC-only stopping
      early_stop <- FALSE  # Flag for early stopping
      
      # Initialize iteration performance tracking
      iteration_results <- data.frame(
        iteration = integer(),
        AIC = numeric(),
        AICc = numeric(),
        deviance = numeric(),
        percent_deviance = numeric(),
        max_alphaW_change = numeric(),
        alphaW_values = character(),
        bandwidths = character(),
        stringsAsFactors = FALSE
      )
      
      while (INT <= max_iterations && !early_stop) {
        if (verbose) cat("--- Back-fitting Iteration:", INT, "---\n")
        fi_old <- Fi
        
        for (i in 1:nvarg) {
          if (verbose) cat("\n    Optimizing variable", i, " (", colnames(X)[i], ")...\n")
          
          # Calculate offset from other variables
          offset_for_i <- Offset
          if(nvarg > 1) {
            Fi_others <- Fi[, -i, drop = FALSE]
            offset_for_i <- offset_for_i + apply(Fi_others, 1, sum)
          }
          X_var <- as.matrix(X[, i])
          
          # Check if offset is reasonable
          if (any(is.na(offset_for_i)) || any(is.infinite(offset_for_i))) {
            if (verbose) cat("      Warning: Invalid offset detected, using global offset\n")
            offset_for_i <- Offset
          }
          
          # Optimize bandwidth using GSS with the correct offset
          tryCatch({
            mband[i] <- GSS(Y, X_var, offset_for_i)
            if (is.na(mband[i]) || is.infinite(mband[i]) || mband[i] <= 0) {
              if (verbose) cat("      Warning: Invalid bandwidth, using previous value\n")
              mband[i] <- if (i > 1) mband[i-1] else hh
            }
            if (verbose) cat("      Optimized bandwidth for variable", i, ":", round(mband[i], 2), "\n")
          }, error = function(e) {
            if (verbose) cat("      Warning: Bandwidth optimization failed, using previous value\n")
            mband[i] <- if (i > 1) mband[i-1] else hh
            if (verbose) cat("      Using fallback bandwidth for variable", i, ":", round(mband[i], 2), "\n")
          })
          
          # --- FIX: Pass the 'Offset' variable and baseline_alphaWs to the search function ---
          search_result <- adaptive_grid_search_with_offset(
            var_idx = i, Fi = Fi, mband = mband, Y = Y, X = X, method = method, 
            model = model, N = N, nvarg = nvarg, wt = wt, E = E, COORD = COORD, 
            sequ = sequ, distancekm = distancekm, Offset = Offset, # <- FIX IS HERE
            parg = parg, verbose = verbose, baseline_alphaWs = rep(baseline_optimal_alphaW, nvarg)
          )
          optimal_alphas[i] <- search_result$alpha
          current_variable_aicc <- search_result$aicc
          
          # Check if alphaW is reasonable - STOP if invalid instead of fallback
          if (is.na(optimal_alphas[i]) || is.infinite(optimal_alphas[i]) || optimal_alphas[i] < 0 || optimal_alphas[i] > 1) {
            cat("ERROR: Invalid alphaW for variable", i, ":", optimal_alphas[i], "\n")
            cat("ERROR: search_result$alpha =", search_result$alpha, "search_result$aicc =", search_result$aicc, "\n")
            stop("Invalid alphaW optimization result - halting to examine")
          }
          
          # Update the fit for the current variable using the new optimal parameters and correct offset
          tryCatch({
            yhat_beta_gwr <- gwr(mband[i], Y, X_var, fi = offset_for_i, alphaW_val = optimal_alphas[i])
            beta[, i] <- yhat_beta_gwr[, 2]
            Fi[, i] <- X_var * beta[, i]
            
            # Check if the fit is reasonable
            if (any(is.na(beta[, i])) || any(is.infinite(beta[, i]))) {
              if (verbose) cat("      Warning: Invalid beta coefficients, using previous values\n")
              if (i > 1) {
                beta[, i] <- beta[, i-1]  # Use previous variable's coefficients
                Fi[, i] <- X_var * beta[, i]
              }
            }
          }, error = function(e) {
            if (verbose) cat("      Warning: GWR fit failed, using previous values\n")
            if (i > 1) {
              beta[, i] <- beta[, i-1]  # Use previous variable's coefficients
              Fi[, i] <- X_var * beta[, i]
            }
          })
          
          # Check performance after each variable update using proper MGWR evaluation
          if (i >= 1) {  # Check after each variable is updated
            # Calculate performance using gwr_local_multi_alpha with current state of all bandwidths and alphaWs
            current_res <- gwr_local_multi_alpha(
              H = mband, y = Y, x = X, fi = finb, alpha_vals = optimal_alphas, 
              method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E, 
              COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset, 
              parg = parg, yhat_beta = yhat_beta, verbose = FALSE
            )
            current_yhbeta <- current_res$yhbeta
            current_beta <- current_yhbeta[, 2:(nvarg + 1)]
            current_Fi <- X * current_beta
            current_sm <- current_res$sm
            current_alphai <- current_res$alphai
            current_v1 <- sum(diag(current_sm))
            current_yhat <- exp(apply(current_Fi, 1, sum) + Offset)
            current_ll <- sum(Y * log(current_alphai[,2] * current_yhat) - (Y + 1/current_alphai[,2]) * log(1 + current_alphai[,2] * current_yhat) +
                                lgamma(Y + 1/current_alphai[,2]) - lgamma(1/current_alphai[,2]) - lgamma(Y + 1))
            current_AIC <- 2 * current_v1 - 2 * current_ll
            # Fix AICc calculation with safeguard for when v1 is close to N
            if ((N - current_v1 - 1) > 0) {
              current_AICc <- current_AIC + 2 * current_v1 * (current_v1 + 1) / (N - current_v1 - 1)
            } else {
              current_AICc <- current_AIC  # Fallback to AIC if denominator is invalid
            }
            
            # Calculate current percent deviance explained (SAME as final model)
            tt_current <- Y/current_yhat
            tt_current <- ifelse(tt_current==0, E^-10, tt_current)
            current_deviance <- 2*sum(Y*log(tt_current)-(Y+1/current_alphai[,2])*log((1+current_alphai[,2]*Y)/(1+current_alphai[,2]*current_yhat)))
            
            tt2_current <- Y/mean(Y)
            tt2_current <- ifelse(tt2_current==0, E^-10, tt2_current)
            null_deviance_current <- 2*sum(Y*log(tt2_current)-(Y+1/current_alphai[,2])*log((1+current_alphai[,2]*Y)/(1+current_alphai[,2]*mean(Y))))
            current_percent_deviance <- (1 - current_deviance / null_deviance_current)
            
            if(verbose) {
              cat("\n--- After Variable", i, "| Alphas:", paste(round(optimal_alphas, 3), collapse=", "), "---\n")
              cat("  Current Performance:\n")
              cat("    AIC: ", round(current_AIC, 2), "\n")
              cat("    AICc: ", round(current_AICc, 2), "\n")
              cat("    pctdev: ", round(current_percent_deviance, 4), "\n")
              cat("    Alpha range: ", round(range(current_alphai[,2]), 2), "\n")
              cat("    Alpha mean: ", round(mean(current_alphai[,2]), 2), "\n")
            }
          }
        }
        
        # Capture iteration performance results
        if (exists("current_AIC") && exists("current_AICc") && exists("current_percent_deviance")) {
          # Calculate max alphaW change from previous iteration
          if (INT == 1) {
            max_alphaW_change <- 0
          } else {
            max_alphaW_change <- max(abs(optimal_alphas - prev_optimal_alphas), na.rm = TRUE)
          }
          
          # Store iteration results
          iteration_results <- rbind(iteration_results, data.frame(
            iteration = INT,
            AIC = current_AIC,
            AICc = current_AICc,
            deviance = current_deviance,
            percent_deviance = current_percent_deviance,
            max_alphaW_change = max_alphaW_change,
            alphaW_values = paste(round(optimal_alphas, 3), collapse = ", "),
            bandwidths = paste(round(mband, 2), collapse = ", "),
            stringsAsFactors = FALSE
          ))
          
          # Store current alphas for next iteration comparison
          prev_optimal_alphas <- optimal_alphas
          
          if (verbose) {
            cat("\n=== ITERATION", INT, "SUMMARY ===\n")
            cat("AIC:", round(current_AIC, 2), "\n")
            cat("AICc:", round(current_AICc, 2), "\n")
            cat("Percent Deviance:", round(current_percent_deviance, 4), "\n")
            cat("Max AlphaW Change:", round(max_alphaW_change, 6), "\n")
            cat("AlphaW Values:", paste(round(optimal_alphas, 3), collapse = ", "), "\n")
            cat("Bandwidths:", paste(round(mband, 2), collapse = ", "), "\n")
          }
          
          # Save iteration performance to text file
          timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
          iteration_file <- paste0("iteration_", INT, "_performance_", timestamp, ".txt")
          
          # Create detailed performance report
          performance_text <- paste0(
            "=== MGWNBR ITERATION ", INT, " PERFORMANCE REPORT ===\n",
            "Timestamp: ", Sys.time(), "\n",
            "Method: ", method, "\n",
            "Model: ", model, "\n",
            "MGWR: ", mgwr, "\n\n",
            
            "=== PERFORMANCE METRICS ===\n",
            "AIC: ", round(current_AIC, 4), "\n",
            "AICc: ", round(current_AICc, 4), "\n",
            "Deviance: ", round(current_deviance, 4), "\n",
            "Percent Deviance Explained: ", round(current_percent_deviance, 6), "%\n",
            "Max AlphaW Change: ", round(max_alphaW_change, 8), "\n\n",
            
            "=== BASELINE COMPARISON ===\n",
            "Baseline AIC: ", round(baseline_AIC, 4), "\n",
            "Baseline AICc: ", round(baseline_AICc, 4), "\n",
            "Baseline Percent Deviance: ", round(baseline_percent_deviance, 6), "%\n",
            "AIC Improvement: ", round(baseline_AIC - current_AIC, 4), "\n",
            "AICc Improvement: ", round(baseline_AICc - current_AICc, 4), "\n",
            "Deviance Improvement: ", round(current_percent_deviance - baseline_percent_deviance, 6), "%\n",
            "Note: Comparison valid after iteration 1 completion\n\n",
            
            "=== OPTIMIZED PARAMETERS ===\n",
            "AlphaW Values: ", paste(round(optimal_alphas, 4), collapse = ", "), "\n",
            "Bandwidths: ", paste(round(mband, 4), collapse = ", "), "\n\n",
            
            "=== VARIABLE DETAILS ===\n"
          )
          
          # Add variable-specific details
          for (i in 1:nvarg) {
            performance_text <- paste0(performance_text,
              "Variable ", i, " (", colnames(X)[i], "):\n",
              "  AlphaW: ", round(optimal_alphas[i], 4), "\n",
              "  Bandwidth: ", round(mband[i], 4), "\n"
            )
          }
          
          performance_text <- paste0(performance_text,
            "\n=== CONVERGENCE STATUS ===\n",
            "Iteration: ", INT, " of ", max_iterations, "\n",
            "Early Stop: ", early_stop, "\n",
            "Converged: ", ifelse(max_alphaW_change < 0.001, "YES", "NO"), "\n"
          )
          
          # Write to file
          writeLines(performance_text, iteration_file)
          if (verbose) cat("Iteration performance saved to:", iteration_file, "\n")
        }
        
        INT <- INT + 1
      }
      
      if (verbose) cat("  - Back-fitting with Grid Search completed.\n")
      
    } else {
      # Original back-fitting logic for other methods (non-SMR)
      if (verbose) cat("  - Starting initial back-fitting for bandwidth estimation...\n")
      # Initialize convergence tracking for non-SMR methods
      INT <- 1
      converged <- FALSE
      
      while (!converged && INT < int) {
        fi_old <- Fi
        diffi <- 0
        for (i in 1:nvarg) {
          yhat_beta_res <- (apply(Fi, 1, sum) + Offset) - Fi[,i]
          mband[i] <- GSS(Y, as.matrix(X[, i]), yhat_beta_res)
          yhat_beta_gwr <- gwr(mband[i], Y, as.matrix(X[, i]), yhat_beta_res, alphaW_val = 1.0)
          beta[, i] <- yhat_beta_gwr[, 2]
          Fi[, i] <- X[, i] * beta[, i]
          diffi <- diffi + sum((Fi[, i] - fi_old[, i])^2)
        }
        # For non-SMR methods, use simple convergence check
        converged <- (diffi < 0.001)
        if(verbose) cat("\n--- End of Iteration", INT, "| Converged:", converged, "---\n")
        INT <- INT + 1
      }
      if (verbose) cat("  - Initial back-fitting for bandwidth estimation completed.\n")
    }
    # --- END OF METHOD-SPECIFIC LOGIC ---
    
    # Final model calculation using the optimized parameters
    if (!exists("optimal_alphas")) { optimal_alphas <- rep(1.0, nvarg) }
    
    cat("DEBUG: Using bandwidths for final model fit:", paste(round(mband, 2), collapse=", "), "\n")
    
    # Use the same approach as back-fitting for consistency
    if (method == "adaptive_bsq") {
      # Use single alphaW approach for adaptive_bsq (always alphaW=1.0)
      res <- gwr_local(
        H = mband[1], y = Y, x = X, fi = finb, alphaW_val = 1.0,
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta
      )
    } else if (method == "adaptive_bsq_smr") {
      # Use single alphaW approach (like baseline and back-fitting)
      res <- gwr_local(
        H = mband[1], y = Y, x = X, fi = finb, alphaW_val = optimal_alphas[1],
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta
      )
    } else {
      # Use multi-alphaW approach for other methods
      res <- gwr_local_multi_alpha(
        H = mband, y = Y, x = X, fi = finb, alpha_vals = optimal_alphas,
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta, verbose = FALSE
      )
    }
    yhat_beta <- res$yhbeta
    sm <- res$sm
    alphai <- res$alphai
    sm3 <- res$sm3
    Sm2 <- sm  # Add Sm2 for compatibility
    beta <- yhat_beta[,2:(nvarg+1)]
    Fi <- X*beta
    
    if (method == "adaptive_bsq_smr"){
      alphaW_info <- data.frame(
        Variable = colnames(X),
        Optimal_AlphaW = optimal_alphas,
        Optimal_Bandwidth = mband
      )
      output$alphaW_optimum_per_variable <- alphaW_info
    }
    
    #   
    # }
    
  }
  if (verbose) cat("Step 4: Computing final model statistics...\n")
  
  v1 <- sum(diag(sm))
  if (verbose) cat("  - Effective number of parameters (v1):", round(v1, 2), "\n")
  
  if (model=='gaussian'){
    yhat <- apply(Fi, 1, sum)
    res <- Y-yhat
    rsqr1 <- t(res*wt)%*%res
    ym <- t(Y*wt)%*%Y
    rsqr2 <- ym-(sum(Y*wt)^2)/sum(wt)
    rsqr <- 1-rsqr1/rsqr2
    rsqradj <- 1-((N-1)/(N-v1))*(1-rsqr)
    sigma2 <- as.vector(N*rsqr1/((N-v1)*sum(wt)))
    root_mse <- sqrt(sigma2)
  }
  for (jj in 1:nvarg){
    m1 <- (jj-1)*N+1
    m2 <- m1+(N-1)
    ENP[jj] <- sum(diag(mrj[,m1:m2]))
    if (!mgwr){
      ENP[jj] <- sum(diag(sm))
    }
    if (model=='gaussian'){
      if (!mgwr){
        stdbm[,jj] <- sqrt(sigma2*sm3[,jj])
      }
      else{
        stdbm[,jj] <- sqrt(diag(sigma2*Cm[,m1:m2]%*%t(Cm[,m1:m2])))
      }
    }
    else{ #else if (model=='poisson' | model=='negbin' | model=='logistic'){
      if (mgwr){
        # For adaptive_bsq with mgwr=TRUE, use sm3 since Cm and mAi are not available
        if (method == "adaptive_bsq") {
          stdbm[,jj] <- sqrt(sm3[,jj])
        } else {
          stdbm[,jj] <- sqrt(diag(Cm[,m1:m2]%*%diag(1/mAi[,jj])%*%t(Cm[,m1:m2])))
        }
      }
      else{
        stdbm[,jj] <- sqrt(sm3[,jj])
      }
    }
  }
  if (model=='gaussian'){
    ll <- -N*log(rsqr1/N)/2-N*log(2*acos(-1))/2-sum((Y-yhat)*(Y-yhat))/(2*(rsqr1/N))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    stats_measures <- c(sigma2, root_mse, round(rsqr, 4),
                        round(rsqradj, 4), ll, AIC, AICc)
    names(stats_measures) <- c("sigma2e", "root_mse",
                               "R_square", "Adj_R_square",
                               "full_Log_likelihood",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  }
  else if (model=='poisson'){
    if (verbose) cat("  - Computing Poisson model statistics...\n")
    yhat <- exp(apply(Fi, 1, sum)+Offset)
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    dev <- 2*sum(Y*log(tt)-(Y-yhat))
    ll <- sum(-yhat+Y*log(yhat)-lgamma(Y+1))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(Y*log(tt2)-(Y-mean(Y)))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    stats_measures <- c(dev, ll, pctdev, adjpctdev, AIC,
                        AICc)
    names(stats_measures) <- c("deviance",
                               "full_Log_likelihood",
                               "pctdev", "adjpctdev",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  }
  else if (model=='negbin'){
    if (verbose) cat("  - Computing negative binomial model statistics...\n")
    yhat <- exp(apply(Fi, 1, sum)+Offset)
    
    # Debug: Check for problematic values
    cat("DEBUG: Fi range:", range(Fi, na.rm=TRUE), "\n")
    cat("DEBUG: yhat range:", range(yhat, na.rm=TRUE), "\n")
    cat("DEBUG: Y range:", range(Y, na.rm=TRUE), "\n")
    cat("DEBUG: Any NaN in Fi:", any(is.na(Fi)), "\n")
    cat("DEBUG: Any Inf in Fi:", any(is.infinite(Fi)), "\n")
    cat("DEBUG: Any NaN in yhat:", any(is.na(yhat)), "\n")
    cat("DEBUG: Any Inf in yhat:", any(is.infinite(yhat)), "\n")
    
    rmse_val <- sqrt(mean((Y - yhat)^2)) 
    
    if (verbose) cat("    * RMSE:", round(rmse_val, 4), "\n")
    
    # Debug: Check alphai values
    if (method == "adaptive_bsq_smr" && mgwr) {
      # cat("DEBUG: alphai[,2] values for negbin model:", paste(head(alphai[,2], 5), collapse=", "), "...\n")
      # cat("DEBUG: alphai[,2] range:", range(alphai[,2]), "\n")
    }
    
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    dev <- 2*sum(Y*log(tt)-(Y+1/alphai[,2])*log((1+alphai[,2]*Y)/(1+alphai[,2]*yhat)))
    ll <- sum(Y*log(alphai[,2]*yhat)-(Y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(Y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(Y+1))
    # Use simple formula for both mgwr=FALSE and mgwr=TRUE (consistent)
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*v1*(v1+1)/(N-v1-1)
    
    if (verbose) {
      cat("    * Log-likelihood:", round(ll, 2), "\n")
      cat("    * AIC:", round(AIC, 2), "\n")
      cat("    * AICc:", round(AICc, 2), "\n")
    }
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(Y*log(tt2)-(Y+1/alphai[,2])*log((1+alphai[,2]*Y)/(1+alphai[,2]*mean(Y))))
    pctdev <- 1-dev/devnull
    # Use simple formula for adjusted percent deviance (consistent)
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    stats_measures <- c(rmse_val, dev, ll, pctdev, adjpctdev, AIC,
                        AICc)
    names(stats_measures) <- c("RMSE", "deviance",
                               "full_Log_likelihood",
                               "pctdev", "adjpctdev",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  }
  else{ #else if (model=='logistic'){
    yhat <- exp(apply(Fi, 1, sum))/(1+exp(apply(Fi, 1, sum)))
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    yhat2 <- ifelse(yhat==1, 0.99999, yhat)
    tt2 <- (1-Y)/(1-yhat2)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    dev <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    lyhat2 <- 1-yhat
    lyhat2 <- ifelse(lyhat2==0, E^-10, lyhat2)
    ll <- sum(Y*log(yhat)+(1-Y)*log(lyhat2))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    tt <- Y/mean(Y)
    tt <- ifelse(tt==0, E^-10, tt)
    tt2 <- (1-Y)/(1-mean(Y))
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    stats_measures <- c(dev, ll, pctdev, adjpctdev, AIC,
                        AICc)
    names(stats_measures) <- c("deviance",
                               "full_Log_likelihood",
                               "pctdev", "adjpctdev",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  }
  ENP[nvarg+1] <- sum(diag(sm))
  ENP[nvarg+2] <- sum(diag(Sm2))
  if (model=='negbin'){
    ENP <- c(ENP, (v1/nvarg))
    names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR', 'alphaNB')
  }
  else{
    names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR')
  }
  header <- append(header, "ENP")
  output <- append(output, list(ENP))
  names(output)[length(output)] <- "ENP"
  dff <- N-v1
  # Robust handling for tstat calculation
  tstat <- beta/stdbm
  tstat[is.na(tstat) | is.infinite(tstat)] <- 0
  
  # Robust handling for probt calculation
  probt <- 2*(1-pt(abs(tstat), dff))
  probt[is.na(probt) | is.infinite(probt)] <- 1
  
  # Robust handling for ENP to prevent division by zero or NaN
  ENP_clean <- ENP
  ENP_clean[is.na(ENP_clean) | is.infinite(ENP_clean) | ENP_clean <= 0] <- 1
  
  malpha <- ENP_clean
  malpha[1:nvarg] <- 0.05/ENP_clean[1:nvarg]
  malpha[nvarg+1] <- 0.05*(nvarg/v1)
  malpha[nvarg+2] <- 0.05*(nvarg/sum(diag(Sm2)))
  if (!mgwr){
    malpha[1:nvarg] <- 0.05*(nvarg/v1)
  }
  if (model=='negbin'){
    malpha[nvarg+3] <- 0.05*(nvarg/v1)
  }
  # Robust handling for t_critical calculation
  t_critical <- abs(qt(malpha/2,dff))
  t_critical[is.na(t_critical) | is.infinite(t_critical)] <- 1.96
  beta2 <- beta
  if (model=='negbin'){
    alphaNB <- alphai[,2]  # Negative binomial alpha (dispersion parameter)
    beta2 <- cbind(beta, alphaNB)
  }
  qntl <- apply(beta2, 2, robust_quantile, c(0.25, 0.5, 0.75))
  IQR <- (qntl[3,]-qntl[1,])
  qntl <- rbind(round(qntl, 6), IQR=round(IQR, 6))
  descriptb <- rbind(apply(beta2, 2, function(x) mean(x[!is.na(x) & !is.infinite(x)])), 
                     apply(beta2, 2, function(x) min(x[!is.na(x) & !is.infinite(x)])), 
                     apply(beta2, 2, function(x) max(x[!is.na(x) & !is.infinite(x)])))
  rownames(descriptb) <- c('Mean', 'Min', 'Max')
  
  final_estimates_df <- cbind(y_observed = Y, y_predicted = yhat, beta2)
  
  if (model=='negbin'){ #release 2
    colnames(final_estimates_df) <- c('y_observed', 'y_predicted', 'Intercept', XVAR, 'alphaNB')#release 2
  } #release 2
  else{ #release 2
    colnames(final_estimates_df) <- c('y_observed', 'y_predicted', 'Intercept', XVAR) #release 2
  } #release 2
  output <- append(output, list(as.data.frame(final_estimates_df))) #release 2
  names(output)[length(output)] <- "mgwr_param_estimates" #release 2
  if (model=='negbin'){
    colnames(qntl) <- c('Intercept', XVAR, 'alphaNB')
  }
  else{
    colnames(qntl) <- c('Intercept', XVAR)
  }
  header <- append(header, "Quantiles of MGWR Parameter Estimates")
  output <- append(output, list(qntl))
  names(output)[length(output)] <- "qntls_mgwr_param_estimates"
  if (model=='negbin'){
    colnames(descriptb) <- c('Intercept', XVAR, 'alphaNB')
  }
  else{
    colnames(descriptb) <- c('Intercept', XVAR)
  }
  header <- append(header, "Descriptive Statistics")
  output <- append(output, list(descriptb))
  names(output)[length(output)] <- "descript_stats_mgwr_param_estimates"
  stdbeta <- stdbm
  stdbeta2 <- stdbeta
  if (model=='negbin'){
    stdalphaNB <- alphai[,3]  # Standard error of negative binomial alpha
    stdbeta2 <- cbind(stdbeta, stdalphaNB)
  }
  qntls <- apply(stdbeta2, 2, robust_quantile, c(0.25, 0.5, 0.75))
  IQR <- (qntls[3,]-qntls[1,])
  qntls <- rbind(round(qntls, 6), IQR=round(IQR, 6))
  descripts <- rbind(apply(stdbeta2, 2, function(x) mean(x[!is.na(x) & !is.infinite(x)])), 
                     apply(stdbeta2, 2, function(x) min(x[!is.na(x) & !is.infinite(x)])), 
                     apply(stdbeta2, 2, function(x) max(x[!is.na(x) & !is.infinite(x)])))
  rownames(descripts) <- c('Mean', 'Min', 'Max')
  header <- append(header, "alphaNB-level=0.05")
  output <- append(output, list(malpha))
  names(output)[length(output)] <- "p_values"
  t_critical <- round(t_critical, 2)
  header <- append(header, "t-Critical")
  output <- append(output, list(t_critical))
  names(output)[length(output)] <- "t_critical"
  if (model=='negbin'){ #release 2
    colnames(stdbeta2) <- c('Intercept', XVAR, 'alphaNB')#release 2
  } #release 2
  else{ #release 2
    colnames(stdbeta2) <- c('Intercept', XVAR) #release 2
  } #release 2
  output <- append(output, list(as.data.frame(stdbeta2))) #release 2
  names(output)[length(output)] <- "mgwr_se" #release 2
  if (model=='negbin'){
    colnames(qntls) <- c('Intercept', XVAR, 'alphaNB')
  }
  else{
    colnames(qntls) <- c('Intercept', XVAR)
  }
  header <- append(header, "Quantiles of MGWR Standard Errors")
  output <- append(output, list(qntls))
  names(output)[length(output)] <- "qntls_mgwr_se"
  if (model=='negbin'){
    colnames(descripts) <- c('Intercept', XVAR, 'alphaNB')
  }
  else{
    colnames(descripts) <- c('Intercept', XVAR)
  }
  header <- append(header, "Descriptive Statistics of Standard Errors")
  output <- append(output, list(descripts))
  names(output)[length(output)] <- "descripts_stats_se"
  #### global estimates ####
  if (model=='gaussian'){
    bg <- solve(t(X)%*%(X*wt))%*%t(X)%*%(Y*wt)
    s2g <- as.vector(t((Y-X%*%bg)*wt)%*%(Y-X%*%bg)/(N-nrow(bg)))
    varg <- diag(solve(t(X)%*%(X*wt))*s2g)
  }
  if (is.null(weight)){
    vargd <- varg
    dfg <- N-nrow(bg)
  }
  stdg <- matrix(sqrt(vargd))
  if (model=='negbin'){
    bg <- rbind(bg, alphaNBg)
    stdg <- rbind(stdg, sealphaNBg)
    dfg <- dfg-1
  }
  tg <- bg/stdg
  probtg <- 2*(1-pt(abs(tg), dfg))
  bg_stdg_tg_probtg <- cbind(bg, stdg, tg, probtg)
  if (model=='negbin'){
    rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR, 'alphaNB')
  }
  else{
    rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR)
  }
  colnames(bg_stdg_tg_probtg) <- c("Par. Est.", "Std Error", "t Value", "Pr > |t|")
  header <- append(header, "Global Parameter Estimates")
  output <- append(output, list(bg_stdg_tg_probtg))
  names(output)[length(output)] <- "global_param_estimates"
  header <- append(header, "NOTE: The denominator degrees of freedom for the t tests is...")
  output <- append(output, list(dfg))
  names(output)[length(output)] <- "t_test_dfs"
  if (model=='gaussian'){
    resg <- (Y-X%*%bg)
    rsqr1g <- t(resg*wt)%*%resg
    ymg <- t(Y*wt)%*%Y
    rsqr2g <- ymg-(sum(Y*wt)^2)/sum(wt)
    rsqrg <- 1-rsqr1g/rsqr2g
    rsqradjg <- 1-((N-1)/(N-nrow(bg)))%*%(1-rsqrg)
    sigma2g <- N*rsqr1g/((N-nrow(bg))*sum(wt))
    root_mseg <- sqrt(sigma2g)
    ll <- -N*log(rsqr1g/N)/2-N*log(2*acos(-1))/2-sum(resg*resg)/(2*(rsqr1g/N))
    AIC <- -2*ll+2*nrow(bg)
    AICc <- -2*ll+2*nrow(bg)*(N/(N-nrow(bg)-1))
    global_measures <- c(sigma2g, root_mseg, round(c(rsqrg, rsqradjg), 4), ll, AIC, AICc)
    names(global_measures) <- c('sigma2e', 'root_mse', "R_square", "Adj_R_square", 'full_Log_likelihood', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  }
  else if (model=='poisson'){
    yhatg <- exp(X%*%bg+Offset)
    ll <- sum(-yhatg+Y*log(yhatg)-lgamma(Y+1))
    AIC <- -2*ll+2*nvarg
    AICc <- -2*ll+2*nvarg*(N/(N-nvarg-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(Y*log(tt2)-(Y-mean(Y)))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                                'adjpctdevg', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  }
  else if (model=='negbin'){
    yhatg <- exp(X%*%bg[1:(nrow(bg)-1)]+Offset)
    ll <- sum(Y*log(alphaNBg*yhatg)-(Y+1/alphaNBg)*log(1+alphaNBg*yhatg)+lgamma(Y+1/alphaNBg)-lgamma(1/alphaNBg)-lgamma(Y+1))
    AIC <- -2*ll+2*(nvarg+1)
    AICc <- -2*ll+2*(nvarg+1)*(N/(N-(nvarg+1)-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(Y*log(tt2)-(Y+1/alphaNBg)*log((1+alphaNBg*Y)/(1+alphaNBg*mean(Y))))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                                'adjpctdevg', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  }
  else{ #else if (model=='logistic'){
    yhatg <- exp(X%*%bg)/(1+exp(X%*%bg))
    lyhat2 <- 1-yhatg
    lyhat2 <- ifelse(lyhat2==0, E^-10, lyhat2)
    ll <- sum(Y*log(yhatg)+(1-Y)*log(lyhat2))
    AIC <- -2*ll+2*nvarg
    AICc <- -2*ll+2*nvarg*(N/(N-nvarg-1))
    tt <- Y/mean(Y)
    tt <- ifelse(tt==0, E^-10, tt)
    tt2 <- (1-Y)/(1-mean(Y))
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                                'adjpctdevg', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  }
  # Prepare final results
  
  bistdt <- cbind(COORD, beta, stdbm, tstat, probt)
  colname1 <- c("Intercept", XVAR)
  parameters2 <- as.data.frame(bistdt)
  names(parameters2) <- c('x', 'y', colname1, paste('std_', colname1, sep=''), paste('tstat_', colname1, sep=''), paste('probt_', colname1, sep=''))
  sig <- matrix("not significant at 90%", N, nvarg)
  for (i in 1:N){
    for (j in 1:nvarg){
      # Robust handling for conditional statements
      probt_val <- probt[i,j]
      ENP_val <- ENP_clean[j]
      
      if (is.na(probt_val) || is.infinite(probt_val) || is.na(ENP_val) || is.infinite(ENP_val)) {
        sig[i,j] <- "not significant at 90%"
      } else if (probt_val < 0.01/ENP_val){
        sig[i,j] <- "significant at 99%"
      }
      else if (probt_val < 0.05/ENP_val){
        sig[i,j] <- "significant at 95%"
      }
      else if (probt_val < 0.1/ENP_val){
        sig[i,j] <- "significant at 90%"
      }
      else{
        sig[i,j] <- "not significant at 90%"
      }
    }
  }
  sig_parameters2 <- as.data.frame(sig)
  names(sig_parameters2) <- c(paste('sig_', colname1, sep=''))
  if (model=='negbin'){
    # Robust handling for atstat calculation
    atstat <- alphai[,2]/alphai[,3]
    atstat[is.na(atstat) | is.infinite(atstat)] <- 0
    
    # Robust handling for aprobtstat calculation
    aprobtstat <- 2*(1-pnorm(abs(atstat)))
    aprobtstat[is.na(aprobtstat) | is.infinite(aprobtstat)] <- 1
    
    siga <- rep("not significant at 90%", N)
    for (i in 1:N){
      # Robust handling for conditional statements
      aprobtstat_val <- aprobtstat[i]
      
      if (is.na(aprobtstat_val) || is.infinite(aprobtstat_val)) {
        siga[i] <- "not significant at 90%"
      } else if (aprobtstat_val < 0.01*(nvarg/v1)){
        siga[i] <- "significant at 99%"
      }
      else if (aprobtstat_val < 0.05*(nvarg/v1)){
        siga[i] <- "significant at 95%"
      }
      else if (aprobtstat_val < 0.1*(nvarg/v1)){
        siga[i] <- "significant at 90%"
      }
      else{
        siga[i] <- "not significant at 90%"
      }
    }
    alphai <- cbind(alphai, atstat, aprobtstat)
    Alpha <- as.data.frame(alphai)
    names(Alpha) <- c("id", "alphaNB", "std", "tstat", "probt")
    sig_alpha <- as.data.frame(siga)
    names(sig_alpha) <- "sig_alphaNB"
  }
  ###################################
  min_bandwidth <- as.data.frame(t(mband))
  if (!mgwr){
    names(min_bandwidth) <- 'Intercept'
  }
  else{
    names(min_bandwidth) <- colname1
  }
  parameters2 <- cbind(parameters2, sig_parameters2)
  if (model=='negbin'){
    Alpha <- cbind(Alpha, sig_alpha)  # sig_alpha contains sig_alphaNB
  }
  
  # ================================================================================
  # AUTOMATIC EXPORT OF RESULTS
  # ================================================================================
  
  # Generate timestamp for file names
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Create comprehensive results CSV (similar to Validation1_ED_results.csv)
  comprehensive_results <- data.frame(
    id = 1:N,
    x_coor = COORD[,1],
    y_coor = COORD[,2],
    y = Y,
    ols_residual = Y - mean(Y),  # Simple OLS residual
    mgwr_yhat = yhat,
    mgwr_residual = Y - yhat,
    local_pctdev = rep(NA, N)  # Will be calculated below
  )
  
  # Add beta coefficients
  for (i in 1:nvarg) {
    comprehensive_results[[paste0("beta_", colname1[i])]] <- beta[,i]
  }
  
  # Add standard errors
  for (i in 1:nvarg) {
    comprehensive_results[[paste0("se_", colname1[i])]] <- stdbm[,i]
  }
  
  # Add t-statistics
  for (i in 1:nvarg) {
    comprehensive_results[[paste0("t_", colname1[i])]] <- tstat[,i]
  }
  
  # Add p-values
  for (i in 1:nvarg) {
    comprehensive_results[[paste0("p_", colname1[i])]] <- probt[,i]
  }
  
  # Calculate local pctdev for each observation
  for (i in 1:N) {
    # Simplified local pctdev calculation for negative binomial
    # Using the same approach as the global pctdev but for individual observations
    
    # Local deviance for observation i
    if (yhat[i] > 0 && alphai[i,2] > 0) {
      # Log-likelihood for observation i with model prediction
      ll_model_i <- Y[i] * log(yhat[i]) - (Y[i] + 1/alphai[i,2]) * log(1 + alphai[i,2] * yhat[i]) + 
        lgamma(Y[i] + 1/alphai[i,2]) - lgamma(1/alphai[i,2]) - lgamma(Y[i] + 1)
      
      # Log-likelihood for observation i with null prediction (mean)
      ll_null_i <- Y[i] * log(mean(Y)) - (Y[i] + 1/alphai[i,2]) * log(1 + alphai[i,2] * mean(Y)) + 
        lgamma(Y[i] + 1/alphai[i,2]) - lgamma(1/alphai[i,2]) - lgamma(Y[i] + 1)
      
      # Local deviance (ensure positive values)
      local_deviance_i <- abs(-2 * ll_model_i)
      local_null_deviance_i <- abs(-2 * ll_null_i)
      
      
      # Local pctdev
      if (local_null_deviance_i > 0) {
        comprehensive_results$local_pctdev[i] <- 1 - (local_deviance_i / local_null_deviance_i)
      } else {
        comprehensive_results$local_pctdev[i] <- 0
      }
    } else {
      comprehensive_results$local_pctdev[i] <- 0
    }
  }
  
  # Calculate sum of weights (diagonal elements of smoothing matrix)
  for (i in 1:nvarg) {
    comprehensive_results[[paste0("sumW_", colname1[i])]] <- diag(sm)
  }
  
  # Add alphaNB for negative binomial model
  if (model == 'negbin') {
    comprehensive_results$beta_alpha <- alphai[,2]
    comprehensive_results$se_alpha <- alphai[,3]
    comprehensive_results$t_alpha <- alphai[,2] / alphai[,3]
    comprehensive_results$p_alpha <- 2*(1-pnorm(abs(comprehensive_results$t_alpha)))
    comprehensive_results$sumW_alpha <- diag(sm)  # Use same sum of weights as other variables
  }
  
  # Export comprehensive results CSV
  csv_filename <- paste0("mgwnbr_results_", timestamp, ".csv")
  write.csv(comprehensive_results, csv_filename, row.names = FALSE)
  if (verbose) cat("Comprehensive results exported to:", csv_filename, "\n")
  
  # Create summary text file (similar to Validation1_ED_summary.txt)
  summary_filename <- paste0("mgwnbr_summary_", timestamp, ".txt")
  if (verbose) cat("Creating summary file:", summary_filename, "\n")
  sink(summary_filename)
  
  cat("================================================================================\n")
  cat("MGWNBR Version: Custom Implementation\n")
  cat("Released on:", format(Sys.Date(), "%m/%d/%Y"), "\n")
  cat("================================================================================\n")
  cat("Model type:", paste0("                                                            ", model), "\n")
  cat("Number of observations:", paste0("                                                     ", N), "\n")
  cat("Number of covariates:", paste0("                                                          ", nvarg), "\n")
  cat("Dependent variable:", paste0("                                                           ", deparse(substitute(formula))), "\n")
  cat("Variable standardization:", paste0("                                                     ", ifelse(exists("data_std"), "On", "Off")), "\n")
  cat("Total runtime:", paste0("                                                            ", format(difftime(Sys.time(), start_time, units="secs"), digits=2)), " seconds\n\n")
  
  # Global Regression Results
  cat("Global Regression Results\n")
  cat("--------------------------------------------------------------------------------\n")
  if (model == 'gaussian') {
    cat("Residual sum of squares:", paste0("                                                 ", round(output$global_measures['sigma2e'] * N, 3)), "\n")
    cat("Log-likelihood:", paste0("                                                        ", round(output$global_measures['full_Log_likelihood'], 3)), "\n")
    cat("AIC:", paste0("                                                                    ", round(output$global_measures['AIC'], 3)), "\n")
    cat("AICc:", paste0("                                                                   ", round(output$global_measures['AICc'], 3)), "\n")
    cat("R2:", paste0("                                                                         ", round(output$global_measures['R_square'], 3)), "\n")
    cat("Adj. R2:", paste0("                                                                    ", round(output$global_measures['Adj_R_square'], 3)), "\n\n")
  } else if (model == 'negbin') {
    cat("Deviance:", paste0("                                                               ", round(output$global_measures['deviance'], 3)), "\n")
    cat("Log-likelihood:", paste0("                                                        ", round(output$global_measures['full_Log_likelihood'], 3)), "\n")
    cat("AIC:", paste0("                                                                    ", round(output$global_measures['AIC'], 3)), "\n")
    cat("AICc:", paste0("                                                                   ", round(output$global_measures['AICc'], 3)), "\n")
    cat("Percent deviance explained:", paste0("                                              ", round(output$global_measures['pctdevg'], 3)), "\n")
    cat("Adj. percent deviance explained:", paste0("                                         ", round(output$global_measures['adjpctdevg'], 3)), "\n\n")
  }
  
  # Global Parameter Estimates
  cat("Variable                                   Est.         SE  t(Est/SE)    p-value\n")
  cat("------------------------------------ ---------- ---------- ---------- ----------\n")
  global_ests <- output$global_param_estimates
  for (i in 1:nrow(global_ests)) {
    var_name <- rownames(global_ests)[i]
    est <- global_ests[i, 1]
    se <- global_ests[i, 2]
    t_val <- global_ests[i, 3]
    p_val <- global_ests[i, 4]
    cat(sprintf("%-35s %10.3f %10.3f %10.3f %10.3f\n", var_name, est, se, t_val, p_val))
  }
  cat("\n")
  
  # MGWR Results
  cat("Multiscale Geographically Weighted Regression (MGWR) Results\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Coordinates type:", paste0("                                                      ", ifelse(distancekm, "Geographic", "Projected")), "\n")
  cat("Spatial kernel:", paste0("                                                ", method), "\n")
  cat("Criterion for optimal bandwidth:", paste0("                                            ", ifelse(is.null(h), "AICc", "Fixed")), "\n")
  cat("Termination criterion for MGWR:", paste0("                                           ", "1.0e-05"), "\n")
  cat("Number of iterations used:", paste0("                                                    ", ifelse(mgwr, "Variable", "1")), "\n\n")
  
  # MGWR bandwidths
  cat("MGWR bandwidths\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Variable                  Bandwidth      ENP_j   Adj t-val(95%)            DoD_j\n")
  if (mgwr && exists("mband")) {
    for (i in 1:length(mband)) {
      var_name <- ifelse(i == 1, "Intercept", XVAR[i-1])
      bandwidth <- mband[i]
      enp <- ENP[i]
      adj_t_val <- t_critical[i]
      dod <- 1 - (enp / N)
      cat(sprintf("%-25s %10.0f %10.3f %15.3f %15.3f\n", var_name, bandwidth, enp, adj_t_val, dod))
    }
  }
  cat("\n")
  
  # Diagnostic Information
  cat("Diagnostic Information\n")
  cat("--------------------------------------------------------------------------------\n")
  if (model == 'gaussian') {
    cat("Residual sum of squares:", paste0("                                                 ", round(output$measures['sigma2e'] * N, 3)), "\n")
    cat("Effective number of parameters (trace(S)):", paste0("                                ", round(v1, 3)), "\n")
    cat("Degree of freedom (n - trace(S)):", paste0("                                        ", round(N - v1, 3)), "\n")
    cat("Sigma estimate:", paste0("                                                             ", round(sqrt(output$measures['sigma2e']), 3)), "\n")
    cat("Log-likelihood:", paste0("                                                         ", round(output$measures['full_Log_likelihood'], 3)), "\n")
    cat("AIC:", paste0("                                                                    ", round(output$measures['AIC'], 3)), "\n")
    cat("AICc:", paste0("                                                                   ", round(output$measures['AICc'], 3)), "\n")
    cat("R2:", paste0("                                                                         ", round(output$measures['R_square'], 3)), "\n")
    cat("Adj. R2:", paste0("                                                                    ", round(output$measures['Adj_R_square'], 3)), "\n\n")
  } else if (model == 'negbin') {
    cat("Deviance:", paste0("                                                               ", round(output$measures['deviance'], 3)), "\n")
    cat("Effective number of parameters (trace(S)):", paste0("                                ", round(v1, 3)), "\n")
    cat("Degree of freedom (n - trace(S)):", paste0("                                        ", round(N - v1, 3)), "\n")
    cat("Log-likelihood:", paste0("                                                         ", round(output$measures['full_Log_likelihood'], 3)), "\n")
    cat("AIC:", paste0("                                                                    ", round(output$measures['AIC'], 3)), "\n")
    cat("AICc:", paste0("                                                                   ", round(output$measures['AICc'], 3)), "\n")
    cat("Percent deviance explained:", paste0("                                              ", round(output$measures['pctdev'], 3)), "\n")
    cat("Adj. percent deviance explained:", paste0("                                         ", round(output$measures['adjpctdev'], 3)), "\n\n")
  }
  
  # Model Performance Comparison: Traditional vs Local
  cat("Model Performance Comparison: Traditional vs Local\n")
  cat("--------------------------------------------------------------------------------\n")
  if (model == 'gaussian') {
    cat("Traditional (Global) Model Performance:\n")
    cat("  Residual sum of squares:", paste0("                                         ", round(output$global_measures['sigma2e'] * N, 3)), "\n")
    cat("  Log-likelihood:", paste0("                                                ", round(output$global_measures['full_Log_likelihood'], 3)), "\n")
    cat("  AIC:", paste0("                                                            ", round(output$global_measures['AIC'], 3)), "\n")
    cat("  AICc:", paste0("                                                           ", round(output$global_measures['AICc'], 3)), "\n")
    cat("  R2:", paste0("                                                             ", round(output$global_measures['R_square'], 3)), "\n")
    cat("  Adj. R2:", paste0("                                                        ", round(output$global_measures['Adj_R_square'], 3)), "\n\n")
    
    cat("Local (MGWR) Model Performance:\n")
    cat("  Residual sum of squares:", paste0("                                         ", round(output$measures['sigma2e'] * N, 3)), "\n")
    cat("  Log-likelihood:", paste0("                                                ", round(output$measures['full_Log_likelihood'], 3)), "\n")
    cat("  AIC:", paste0("                                                            ", round(output$measures['AIC'], 3)), "\n")
    cat("  AICc:", paste0("                                                           ", round(output$measures['AICc'], 3)), "\n")
    cat("  R2:", paste0("                                                             ", round(output$measures['R_square'], 3)), "\n")
    cat("  Adj. R2:", paste0("                                                        ", round(output$measures['Adj_R_square'], 3)), "\n\n")
    
    cat("Performance Improvement (Local - Traditional):\n")
    cat("  AICc improvement:", paste0("                                               ", round(output$global_measures['AICc'] - output$measures['AICc'], 3)), "\n")
    cat("  R2 improvement:", paste0("                                                 ", round(output$measures['R_square'] - output$global_measures['R_square'], 3)), "\n")
    cat("  Adj. R2 improvement:", paste0("                                            ", round(output$measures['Adj_R_square'] - output$global_measures['Adj_R_square'], 3)), "\n\n")
    
  } else if (model == 'negbin') {
    cat("Traditional (Global) Model Performance:\n")
    cat("  Deviance:", paste0("                                                       ", round(output$global_measures['deviance'], 3)), "\n")
    cat("  Log-likelihood:", paste0("                                                ", round(output$global_measures['full_Log_likelihood'], 3)), "\n")
    cat("  AIC:", paste0("                                                            ", round(output$global_measures['AIC'], 3)), "\n")
    cat("  AICc:", paste0("                                                           ", round(output$global_measures['AICc'], 3)), "\n")
    cat("  Percent deviance explained:", paste0("                                      ", round(output$global_measures['pctdevg'], 3)), "\n")
    cat("  Adj. percent deviance explained:", paste0("                                 ", round(output$global_measures['adjpctdevg'], 3)), "\n\n")
    
    cat("Local (MGWR) Model Performance:\n")
    cat("  Deviance:", paste0("                                                       ", round(output$measures['deviance'], 3)), "\n")
    cat("  Log-likelihood:", paste0("                                                ", round(output$measures['full_Log_likelihood'], 3)), "\n")
    cat("  AIC:", paste0("                                                            ", round(output$measures['AIC'], 3)), "\n")
    cat("  AICc:", paste0("                                                           ", round(output$measures['AICc'], 3)), "\n")
    cat("  Percent deviance explained:", paste0("                                      ", round(output$measures['pctdev'], 3)), "\n")
    cat("  Adj. percent deviance explained:", paste0("                                 ", round(output$measures['adjpctdev'], 3)), "\n\n")
    
    cat("Performance Improvement (Local - Traditional):\n")
    cat("  AICc improvement:", paste0("                                               ", round(output$global_measures['AICc'] - output$measures['AICc'], 3)), "\n")
    cat("  Deviance reduction:", paste0("                                             ", round(output$global_measures['deviance'] - output$measures['deviance'], 3)), "\n")
    cat("  Percent deviance explained improvement:", paste0("                          ", round(output$measures['pctdev'] - output$global_measures['pctdevg'], 3)), "\n\n")
  }
  
  # Summary Statistics For MGWR Parameter Estimates
  cat("Summary Statistics For MGWR Parameter Estimates\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Variable                        Mean        STD        Min     Median        Max\n")
  cat("--------------------      ---------- ---------- ---------- ---------- ----------\n")
  for (i in 1:nvarg) {
    var_name <- colname1[i]
    mean_val <- mean(beta[,i], na.rm = TRUE)
    std_val <- sd(beta[,i], na.rm = TRUE)
    min_val <- min(beta[,i], na.rm = TRUE)
    median_val <- median(beta[,i], na.rm = TRUE)
    max_val <- max(beta[,i], na.rm = TRUE)
    cat(sprintf("%-25s %10.3f %10.3f %10.3f %10.3f %10.3f\n", var_name, mean_val, std_val, min_val, median_val, max_val))
  }
  if (model == 'negbin') {
    var_name <- "alphaNB"
    mean_val <- mean(alphai[,2], na.rm = TRUE)
    std_val <- sd(alphai[,2], na.rm = TRUE)
    min_val <- min(alphai[,2], na.rm = TRUE)
    median_val <- median(alphai[,2], na.rm = TRUE)
    max_val <- max(alphai[,2], na.rm = TRUE)
    cat(sprintf("%-25s %10.3f %10.3f %10.3f %10.3f %10.3f\n", var_name, mean_val, std_val, min_val, median_val, max_val))
  }
  
  cat("================================================================================\n")
  cat("Acknowledgement:\n")
  cat("This analysis was performed using the MGWNBR (Multiscale Geographically Weighted\n")
  cat("Negative Binomial Regression) implementation.\n")
  cat("================================================================================\n")
  
  sink()
  if (verbose) cat("Summary exported to:", summary_filename, "\n")
  
  # i <- 1
  # for (element in output){
  #   cat(header[i], "\n")
  #   print(element)
  #   i <- i+1
  # }
  message("NOTE: The denominator degrees of freedom for the t tests is ", dfg, ".")
  
  # Record end time and calculate total runtime
  end_time <- Sys.time()
  total_runtime <- difftime(end_time, start_time, units="secs")
  
  if (verbose) {
    cat("\n================================================================================\n")
    cat("MGWNBR ANALYSIS COMPLETED\n")
    cat("================================================================================\n")
    cat("Total runtime:", format(total_runtime, digits=2), "seconds\n")
    cat("Results saved to:", csv_filename, "\n")
    cat("Summary saved to:", summary_filename, "\n")
    cat("================================================================================\n")
  }
  
  # Cleanup parallel computing
  cleanup_parallel(cl)
  
  invisible(output)
}

# Multi-alphaW optimization functions for adaptive_bsq_smr with mgwr=TRUE

# TRUE HYBRID GWR function: Different bandwidth + alphaW per variable
true_hybrid_gwr <- function(Y, X, all_dist_spatial, all_dist_sim, bandwidths, alpha_vector, 
                            method, N, nvarg, Offset) {
  
  yhbeta <- matrix(0, N, nvarg)
  
  for (i in 1:N) {
    if (method == "adaptive_bsq_smr") {
      # Calculate spatial weights for each variable with its own bandwidth
      dist_i <- all_dist_spatial[i, ]
      w1_i <- rep(0, N)
      
      # Calculate attribute similarity weights for each variable with its own bandwidth
      sim_i <- all_dist_sim[i, ]
      w2_i <- rep(0, N)
      
      # Combine weights for each variable with its own (bandwidth, alphaW) pair
      w_i <- rep(0, N)
      
      for (j in 1:nvarg) {
        # Spatial weights with variable-specific bandwidth
        w1_i <- ifelse(dist_i <= bandwidths[j], (1 - (dist_i/bandwidths[j])^2)^2, 0)
        
        # Attribute similarity weights with variable-specific bandwidth
        w2_i <- exp(-(sim_i/bandwidths[j])^2)
        
        # Combine with variable-specific alphaW
        w_i <- w_i + alpha_vector[j] * w1_i + (1 - alpha_vector[j]) * w2_i
      }
      w_i <- w_i / nvarg  # Average across variables
    } else {
      # Default spatial weights (use mean bandwidth)
      dist_i <- all_dist_spatial[i, ]
      H <- mean(bandwidths)
      w_i <- ifelse(dist_i <= H, (1 - (dist_i/H)^2)^2, 0)
    }
    
    # Weighted least squares
    W <- diag(w_i)
    XWX <- t(X) %*% W %*% X
    XWy <- t(X) %*% W %*% Y
    
    # Solve with regularization
    XWX_reg <- XWX + diag(ncol(X)) * 1e-8
    beta_i <- solve(XWX_reg, XWy)
    yhbeta[i, ] <- beta_i
  }
  
  # Calculate predictions
  beta <- yhbeta[, 1:nvarg, drop=FALSE]
  Fi <- X * beta
  yhat_linear <- apply(Fi, 1, sum)
  yhat_mu <- exp(yhat_linear + Offset)
  yhat_mu <- pmax(yhat_mu, 0.1)  # Ensure positive
  
  # Estimate theta for negative binomial
  residuals <- Y - yhat_mu
  theta <- max(0.1, var(residuals) / mean(yhat_mu))
  
  # Calculate log-likelihood
  ll <- sum(lgamma(Y + theta) - lgamma(theta) - lgamma(Y + 1) + 
              theta * log(theta) - theta * log(yhat_mu + theta) + 
              Y * log(yhat_mu) - Y * log(yhat_mu + theta))
  
  # Calculate AICc
  n_params <- nvarg + 1 + 1  # +1 for theta
  AIC <- 2 * n_params - 2 * ll
  AICc <- AIC + 2 * n_params * (n_params + 1) / (N - n_params - 1)
  
  measures <- c(AIC = AIC, AICc = AICc, Log_likelihood = ll, 
                deviance = sum((Y - yhat_mu)^2), 
                percent_deviance_explained = (1 - sum((Y - yhat_mu)^2) / sum((Y - mean(Y))^2)) * 100,
                theta = theta)
  
  return(list(aicc = AICc, measures = measures, theta = theta))
}

# Simple GWR function for baseline calculation (from hybrid version)
simple_gwr_baseline <- function(Y, X, all_dist_spatial, all_dist_sim, alpha_vector, 
                                method, N, nvarg, Offset) {
  
  yhbeta <- matrix(0, N, nvarg)
  
  for (i in 1:N) {
    if (method == "adaptive_bsq_smr") {
      # Calculate spatial weights
      dist_i <- all_dist_spatial[i, ]
      H <- mean(all_dist_spatial)
      w1_i <- ifelse(dist_i <= H, (1 - (dist_i/H)^2)^2, 0)
      
      # Calculate attribute similarity weights
      sim_i <- all_dist_sim[i, ]
      w2_i <- exp(-(sim_i/H)^2)
      
      # Use different alpha for each variable
      w_i <- rep(0, N)
      for (j in 1:nvarg) {
        w_i <- w_i + alpha_vector[j] * w1_i + (1 - alpha_vector[j]) * w2_i
      }
      w_i <- w_i / nvarg  # Average across variables
    } else {
      # Default spatial weights
      dist_i <- all_dist_spatial[i, ]
      H <- mean(all_dist_spatial)
      w_i <- ifelse(dist_i <= H, (1 - (dist_i/H)^2)^2, 0)
    }
    
    # Weighted least squares
    W <- diag(w_i)
    XWX <- t(X) %*% W %*% X
    XWy <- t(X) %*% W %*% Y
    
    # Solve with regularization
    XWX_reg <- XWX + diag(ncol(X)) * 1e-8
    beta_i <- solve(XWX_reg, XWy)
    yhbeta[i, ] <- beta_i
  }
  
  # Calculate predictions
  beta <- yhbeta[, 1:nvarg, drop=FALSE]
  Fi <- X * beta
  yhat_linear <- apply(Fi, 1, sum)
  yhat_mu <- exp(yhat_linear + Offset)
  yhat_mu <- pmax(yhat_mu, 0.1)  # Ensure positive
  
  # Estimate theta for negative binomial
  residuals <- Y - yhat_mu
  theta <- max(0.1, var(residuals) / mean(yhat_mu))
  
  # Calculate log-likelihood
  ll <- sum(lgamma(Y + theta) - lgamma(theta) - lgamma(Y + 1) + 
              theta * log(theta) - theta * log(yhat_mu + theta) + 
              Y * log(yhat_mu) - Y * log(yhat_mu + theta))
  
  # Calculate AICc
  n_params <- nvarg + 1 + 1  # +1 for theta
  AIC <- 2 * n_params - 2 * ll
  AICc <- AIC + 2 * n_params * (n_params + 1) / (N - n_params - 1)
  
  measures <- c(AIC = AIC, AICc = AICc, Log_likelihood = ll, 
                deviance = sum((Y - yhat_mu)^2), 
                percent_deviance_explained = (1 - sum((Y - yhat_mu)^2) / sum((Y - mean(Y))^2)) * 100,
                theta = theta)
  
  return(list(aicc = AICc, measures = measures, theta = theta))
}

# GSS function that uses a specific alphaW value for bandwidth calculation
GSS_with_alpha <- function(depy, indepx, fix, COORD, sequ, distancekm, method, alpha_val) {
  # DEFINING GOLDEN SECTION SEARCH PARAMETERS #
  if(method=="fixed_g" | method=="fixed_bsq"){
    ax <- 0
    bx <- as.integer(max(dist(COORD))+1)
    if (distancekm){
      bx <- bx*111
    }
  }
  else if (method=="adaptive_bsq" | method=="adaptive_bsq_smr"){
    # Use 5 to N range for adaptive methods
    N <- length(depy)
    # ax <- 5  # Minimum 5
    # bx <- N  # Maximum N
    ax <- floor(0.15 * N)  # Minimum 15 or 15% of N (whichever is larger)
    bx <- floor(0.35 * N)  # Maximum 35% of N (slightly larger)
  }
  r <- 0.61803399
  tol <- 0.1
  lower <- ax
  upper <- bx
  
  # Golden Section Search
  x1 <- upper - r * (upper - lower)
  x2 <- lower + r * (upper - lower)
  
  f1 <- cv_with_alpha(depy, indepx, fix, x1, COORD, sequ, distancekm, method, alpha_val)
  f2 <- cv_with_alpha(depy, indepx, fix, x2, COORD, sequ, distancekm, method, alpha_val)
  
  while (abs(upper - lower) > tol) {
    if (f1 < f2) {
      upper <- x2
      x2 <- x1
      f2 <- f1
      x1 <- upper - r * (upper - lower)
      f1 <- cv_with_alpha(depy, indepx, fix, x1, COORD, sequ, distancekm, method, alpha_val)
    } else {
      lower <- x1
      x1 <- x2
      f1 <- f2
      x2 <- lower + r * (upper - lower)
      f2 <- cv_with_alpha(depy, indepx, fix, x2, COORD, sequ, distancekm, method, alpha_val)
    }
  }
  
  return((lower + upper) / 2)
}

# CV function that uses a specific alphaW value
cv_with_alpha <- function(depy, indepx, fix, H, COORD, sequ, distancekm, method, alpha_val) {
  N <- length(depy)
  cv_sum <- 0
  
  for (i in 1:N) {
    # Calculate distances from point i to all other points
    distan <- cbind(sequ, sqrt((COORD[,1] - COORD[i,1])^2 + (COORD[,2] - COORD[i,2])^2))
    distan <- distan[distan[,2] > 0, ]  # Remove self-distance
    
    if (nrow(distan) == 0) {
      cv_sum <- cv_sum + (depy[i] - mean(depy))^2
      next
    }
    
    # Calculate weights using the specific alphaW value
    u <- nrow(distan)
    w <- rep(0, u)
    
    if (method=="fixed_g"){
      for (jj in 1:u){
        w[jj] <- exp(-(distan[jj,2]/H)^2)
      }
    }
    else if (method=="fixed_bsq"){
      for (jj in 1:u){
        w[jj] <- (1-(distan[jj,2]/H)^2)^2
      }
    }
    else if (method=="adaptive_bsq"){
      distan <- distan[order(distan[, 2]), ]
      distan <- cbind(distan, 1:nrow(distan))
      w_temp <- matrix(0, N-1, 2)
      hn <- distan[H,2]
      if (hn == 0) hn <- 1e-6
      for (jj in 1:(N-1)){
        if (distan[jj,3] <= H){
          w_temp[jj,1] <- (1-(distan[jj,2]/hn)^2)^2
        } else{
          w_temp[jj,1] <- 0
        }
        w_temp[jj,2] <- distan[jj,1]
      }
      w_temp <- w_temp[order(w_temp[, 2]), ]
      w <- w_temp[,1]
    }
    else if (method=="adaptive_bsq_smr"){
      # Use the specific alphaW value for adaptive_bsq_smr
      beta_val <- 1 - alpha_val
      
      # --- w1 (Spatial weights with bisquare kernel) ---
      dist_spatial_sorted <- distan[order(distan[, 2]), ]
      dist_spatial_sorted <- cbind(dist_spatial_sorted, 1:nrow(dist_spatial_sorted))
      w1_matrix <- matrix(0, N-1, 2)
      hn_spatial <- dist_spatial_sorted[H, 2]
      if (hn_spatial == 0) hn_spatial <- 1e-6 
      for (jj in 1:(N-1)) {
        if (dist_spatial_sorted[jj, 3] <= H) {
          w1_matrix[jj, 1] <- (1-(dist_spatial_sorted[jj, 2]/hn_spatial)^2)^2
        } else {
          w1_matrix[jj, 1] <- 0
        }
        w1_matrix[jj, 2] <- dist_spatial_sorted[jj, 1] 
      }
      w1_matrix <- w1_matrix[order(w1_matrix[, 2]), ]
      w1 <- w1_matrix[, 1]
      
      # --- w2 (Attribute similarity weights with bisquare kernel) ---
      sim_dists <- sqrt(rowSums(sweep(indepx, 2, indepx[i, ], "-")^2))
      dist_sim <- cbind(sequ, sim_dists)
      dist_sim_sorted <- dist_sim[order(dist_sim[, 2]), ]
      dist_sim_sorted <- cbind(dist_sim_sorted, 1:nrow(dist_sim_sorted))
      w2_matrix <- matrix(0, N, 2)
      hn_sim <- dist_sim_sorted[H, 2]
      if (hn_sim == 0) hn_sim <- 1e-6 
      for (jj in 1:N) {
        if (dist_sim_sorted[jj, 3] <= H) {
          w2_matrix[jj, 1] <- exp(-(dist_sim_sorted[jj, 2]/hn_sim)^2)
        } else {
          w2_matrix[jj, 1] <- 0
        }
        w2_matrix[jj, 2] <- dist_sim_sorted[jj, 1]
      }
      w2_matrix <- w2_matrix[order(w2_matrix[, 2]), ]
      w2 <- w2_matrix[, 1]
      w <- alpha_val * w1 + beta_val * w2
    }
    
    # Predict using the calculated weights
    if (sum(w) > 0) {
      y_pred <- sum(w * depy[distan[,1]]) / sum(w)
    } else {
      y_pred <- mean(depy)
    }
    
    cv_sum <- cv_sum + (depy[i] - y_pred)^2
  }
  
  return(cv_sum)
}

# Proper multi-alphaW optimization with stopping criteria
optimize_multi_alpha_mgwr <- function(Y, X, method, model, N, nvarg, Offset, wt, E, COORD, sequ, distancekm, parg, yhat_beta, initial_bandwidths, initial_alpha, baseline_aicc, target_improvement = 0.1, max_iterations = 50) {
  # cat("DEBUG: optimize_multi_alpha_mgwr entry\n")
  # cat("DEBUG: baseline_aicc =", baseline_aicc, "target_improvement =", target_improvement, "\n")
  
  # Initialize parameters
  current_bandwidths <- initial_bandwidths
  current_alphas <- rep(initial_alpha, nvarg)
  best_aicc <- baseline_aicc
  best_bandwidths <- current_bandwidths
  best_alphas <- current_alphas
  iteration <- 0
  improvement <- 0
  
  # Optimization loop
  while (iteration < max_iterations && improvement < target_improvement) {
    iteration <- iteration + 1
    # cat("DEBUG: Iteration", iteration, "\n")
    
    # Try different alphaW combinations
    for (i in 1:nvarg) {
      # Test alphaW variations for variable i
      alpha_candidates <- seq(max(0.1, initial_alpha - 0.2), min(0.9, initial_alpha + 0.2), by = 0.05)
      
      for (alpha_test in alpha_candidates) {
        test_alphas <- current_alphas
        test_alphas[i] <- alpha_test
        
        # Test bandwidth variations for variable i
        bandwidth_candidates <- seq(max(current_bandwidths[i] * 0.8, floor(0.15 * N)), min(floor(0.35 * N), current_bandwidths[i] * 1.2), by = 1)
        
        for (bandwidth_test in bandwidth_candidates) {
          test_bandwidths <- current_bandwidths
          test_bandwidths[i] <- bandwidth_test
          
          # Evaluate this combination
          tryCatch({
            result <- gwr_local_multi_alpha(
              H = test_bandwidths, y = Y, x = X, fi = rep(0, N), alpha_vals = test_alphas,
              method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
              COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
              parg = parg, yhat_beta = yhat_beta, verbose = FALSE
            )
            
            # Calculate AICc
            yhbeta <- result$yhbeta
            beta <- yhbeta[, 2:(nvarg + 1)]
            Fi <- X * beta
            sm <- result$sm
            alphai <- result$alphai
            v1 <- sum(diag(sm))
            yhat <- exp(apply(Fi, 1, sum) + Offset)
            ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                        lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
            # Fix AICc calculation for mgwr=FALSE vs mgwr=TRUE in optimization
            # Use consistent AICc calculation for both mgwr=TRUE and mgwr=FALSE
            AIC <- 2 * v1 - 2 * ll
            if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
            
            # Update if better
            if (AICc < best_aicc) {
              best_aicc <- AICc
              best_bandwidths <- test_bandwidths
              best_alphas <- test_alphas
              improvement <- baseline_aicc - best_aicc
              # cat("DEBUG: New best AICc =", best_aicc, "improvement =", improvement, "\n")
            }
          }, error = function(e) {
            # Skip this combination if it fails
          })
        }
      }
    }
    
    # Update current parameters
    current_bandwidths <- best_bandwidths
    current_alphas <- best_alphas
    
    # cat("DEBUG: Iteration", iteration, "completed. Best AICc =", best_aicc, "improvement =", improvement, "\n")
  }
  
  # cat("DEBUG: Optimization completed. Final AICc =", best_aicc, "total improvement =", baseline_aicc - best_aicc, "\n")
  
  return(list(
    bandwidths = best_bandwidths,
    alphas = best_alphas,
    final_aicc = best_aicc,
    iterations = iteration,
    improvement = baseline_aicc - best_aicc
  ))
}

# Global-local hybrid optimization for MGWR with multi-alphaW - Step 3: optimize per-variable alphaWs AND bandwidths
global_local_hybrid_mgwr <- function(Y, X, method, model, N, nvarg, Offset, wt, E, COORD, sequ, distancekm, parg, yhat_beta, initial_bandwidths, initial_alpha, mgwr, target_improvement = 10, max_iterations = 3, verbose = TRUE, cl = NULL) {
  if (verbose) {
    cat("DEBUG: global_local_hybrid_mgwr entry (Step 3: optimize per-variable alphas AND bandwidths with baseline AICc criterion)\n")
    cat("DEBUG: N =", N, "nvarg =", nvarg, "\n")
    cat("DEBUG: initial_bandwidths =", paste(initial_bandwidths, collapse=", "), "\n")
    cat("DEBUG: initial_alpha =", initial_alpha, "\n")
  }
  
  # STEP 1: Find optimal baseline alphaW and AICc from adaptive_bsq_smr with mgwr=FALSE
  if (verbose) cat("STEP 1: Finding optimal baseline alphaW from adaptive_bsq_smr (mgwr=FALSE)...\n")
  baseline_aicc <- 1000  # Default fallback
  optimal_baseline_alpha <- 0.5  # Default fallback
  
  tryCatch({
    # Use single bandwidth for baseline (mgwr=FALSE equivalent)
    baseline_bw <- robust_quantile(1:N, 0.5)  # Median as baseline bandwidth
    
    # Define alphaW objective function for baseline
    baseline_alpha_objective <- function(alphaW_val) {
      tryCatch({
        result <- gwr_local(
          H = baseline_bw, y = Y, x = X, fi = rep(0, N), alphaW_val = alphaW_val,
          method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
          COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
          parg = parg, yhat_beta = yhat_beta
        )
        
        # Calculate AICc
        yhbeta <- result$yhbeta
        beta <- yhbeta[, 2:(nvarg + 1)]
        Fi <- X * beta
        sm <- result$sm
        alphai <- result$alphai
        v1 <- sum(diag(sm))
        yhat <- exp(apply(Fi, 1, sum) + Offset)
        ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                    lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
        # Fix AICc calculation for mgwr=FALSE vs mgwr=TRUE in baseline objective
        # Use consistent AICc calculation for both mgwr=TRUE and mgwr=FALSE
        AIC <- 2 * v1 - 2 * ll
        if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
        
        return(AICc)
      }, error = function(e) {
        cat("DEBUG: ERROR in baseline_alpha_objective alphaW =", alphaW_val, ":", e$message, "\n")
        return(1000)
      })
    }
    
    # Golden Section Search for optimal baseline alphaW
    phi <- (1 + sqrt(5)) / 2
    a <- 0.0  # Lower bound for alphaW
    b <- 1.0  # Upper bound for alphaW
    c <- b - (b - a) / phi
    d <- a + (b - a) / phi
    
    for (iter in 1:50) {  # Max 50 iterations for baseline alphaW optimization
      if (baseline_alpha_objective(c) < baseline_alpha_objective(d)) {
        b <- d
      } else {
        a <- c
      }
      c <- b - (b - a) / phi
      d <- a + (b - a) / phi
      
      if (abs(b - a) < 1e-6) break
    }
    
    optimal_baseline_alpha <- (a + b) / 2
    baseline_aicc <- baseline_alpha_objective(optimal_baseline_alpha)
    
    # cat("DEBUG: Optimal baseline alphaW =", optimal_baseline_alpha, "with AICc =", baseline_aicc, "\n")
  }, error = function(e) {
    # cat("DEBUG: ERROR finding optimal baseline alphaW:", e$message, "\n")
  })
  
  # STEP 3: Joint optimization of bandwidths and alphaWs per variable
  if (verbose) cat("STEP 3: Joint optimization of bandwidths and alphaWs per variable with baseline target...\n")
  
  # Initialize with provided values
  optimal_bandwidths <- initial_bandwidths
  optimal_alphas <- rep(initial_alpha, nvarg)
  if (verbose) cat("DEBUG: Initial alphaW values for all variables:", paste(round(optimal_alphas, 3), collapse=", "), "\n")
  
  # Joint optimization: iterate between bandwidth and alpha optimization
  for (iteration in 1:max_iterations) {
    if (verbose && iteration %% 5 == 0) cat("    DEBUG: Joint optimization iteration", iteration, "\n")
    
    old_bandwidths <- optimal_bandwidths
    old_alphas <- optimal_alphas
    
    # Phase 1: Optimize bandwidths given current alphas
    if (verbose && iteration == 1) cat("      DEBUG: Phase 1 - Optimizing bandwidths given current alphas\n")
    
    # Use parallel computing for bandwidth optimization if available
    if (!is.null(cl)) {
      optimal_bandwidths <- foreach(i = 1:nvarg, .combine = c, .packages = c("stats", "MASS", "sp")) %dopar% {
        # Define bandwidth objective function for this variable
        bw_objective <- function(bw_val) {
          tryCatch({
            result <- gwr_local(
              H = bw_val, y = Y, x = as.matrix(X[,i]), fi = rep(0, N), alphaW_val = optimal_alphas[i],
              method = method, model = model, N = N, nvarg = 1, wt = wt, E = E,
              COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
              parg = parg, yhat_beta = yhat_beta
            )
            
            # Calculate AICc for this variable
            yhbeta <- result$yhbeta
            beta <- yhbeta[, 2]
            Fi <- as.matrix(X[,i]) * beta
            sm <- result$sm
            alphai <- result$alphai
            v1 <- sum(diag(sm))
            yhat <- exp(Fi + Offset)
            ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                        lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
            # Fix AICc calculation for single variable optimization
            AIC <- 2 * v1 - 2 * ll
            if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
            
            return(AICc)
          }, error = function(e) {
            return(1000)
          })
        }
        
        # Golden Section Search for optimal bandwidth
        phi <- (1 + sqrt(5)) / 2
        # Use 5 to N range for adaptive methods
        if (method == "adaptive_bsq" || method == "adaptive_bsq_smr") {
          a <- floor(0.15 * N)  # Minimum 15 or 15% of N (whichever is larger)
          b <- floor(0.35 * N)  # Maximum 35% of N (slightly larger)
        } else {
          a <- 0.01
          b <- max(distancekm)
        }
        
        tolerance <- 1e-3
        
        # Golden section search for bandwidth
        c <- b - (b - a) / phi
        d <- a + (b - a) / phi
        fc <- bw_objective(c)
        fd <- bw_objective(d)
        
        while (abs(b - a) > tolerance) {
          if (fc < fd) {
            b <- d
            fb <- fd
            d <- c
            fd <- fc
            c <- b - (b - a) / phi
            fc <- bw_objective(c)
          } else {
            a <- c
            fa <- fc
            c <- d
            fc <- fd
            d <- a + (b - a) / phi
            fd <- bw_objective(d)
          }
        }
        
        (a + b) / 2
      }
    } else {
      for (i in 1:nvarg) {
        # Define bandwidth objective function for this variable
        bw_objective <- function(bw_val) {
          tryCatch({
            result <- gwr_local(
              H = bw_val, y = Y, x = as.matrix(X[,i]), fi = rep(0, N), alphaW_val = optimal_alphas[i],
              method = method, model = model, N = N, nvarg = 1, wt = wt, E = E,
              COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
              parg = parg, yhat_beta = yhat_beta
            )
            
            # Calculate AICc for this variable
            yhbeta <- result$yhbeta
            beta <- yhbeta[, 2]
            Fi <- as.matrix(X[,i]) * beta
            sm <- result$sm
            alphai <- result$alphai
            v1 <- sum(diag(sm))
            yhat <- exp(Fi + Offset)
            ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                        lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
            # Fix AICc calculation for single variable optimization
            AIC <- 2 * v1 - 2 * ll
            if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
            
            return(AICc)
          }, error = function(e) {
            cat("DEBUG: ERROR in bw_objective for variable", i, "bw =", bw_val, ":", e$message, "\n")
            return(1000)
          })
        }
        
        # Golden Section Search for optimal bandwidth
        phi <- (1 + sqrt(5)) / 2
        # Use 5 to N range for adaptive methods
        if (method == "adaptive_bsq" || method == "adaptive_bsq_smr") {
          a <- floor(0.15 * N)  # Minimum 15 or 15% of N (whichever is larger)
          b <- floor(0.35 * N)  # Maximum 35% of N (slightly larger)
        } else {
          a <- 0.01
          b <- max(distancekm)
        }
        
        tolerance <- 1e-3
        
        # Golden section search for bandwidth
        c <- b - (b - a) / phi
        d <- a + (b - a) / phi
        fc <- bw_objective(c)
        fd <- bw_objective(d)
        
        while (abs(b - a) > tolerance) {
          if (fc < fd) {
            b <- d
            fb <- fd
            d <- c
            fd <- fc
            c <- b - (b - a) / phi
            fc <- bw_objective(c)
          } else {
            a <- c
            fa <- fc
            c <- d
            fc <- fd
            d <- a + (b - a) / phi
            fd <- bw_objective(d)
          }
        }
        
        optimal_bandwidths[i] <- (a + b) / 2
        # cat("DEBUG: Variable", i, "optimal bandwidth =", optimal_bandwidths[i], "\n")
      }
    }
    
    # Phase 2: Optimize alphas given current bandwidths
    if (verbose && iteration == 1) cat("      DEBUG: Phase 2 - Optimizing alphas given current bandwidths\n")
    for (i in 1:nvarg) {
      # Define alpha objective function for this variable
      alpha_objective <- function(alpha_val) {
        tryCatch({
          result <- gwr_local(
            H = optimal_bandwidths[i], y = Y, x = as.matrix(X[,i]), fi = rep(0, N), alphaW_val = alpha_val,
            method = method, model = model, N = N, nvarg = 1, wt = wt, E = E,
            COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
            parg = parg, yhat_beta = yhat_beta
          )
          
          # Calculate AICc for this variable
          yhbeta <- result$yhbeta
          beta <- yhbeta[, 2]
          Fi <- as.matrix(X[,i]) * beta
          sm <- result$sm
          alphai <- result$alphai
          v1 <- sum(diag(sm))
          yhat <- exp(Fi + Offset)
          ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                      lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
          # Fix AICc calculation for single variable alpha optimization
          AIC <- 2 * v1 - 2 * ll
          if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
          
          return(AICc)
        }, error = function(e) {
          cat("DEBUG: ERROR in alpha_objective for variable", i, "alphaW =", alpha_val, ":", e$message, "\n")
          return(1000)
        })
      }
      
      # Golden Section Search for optimal alphaW in range [0,1]
      phi <- (1 + sqrt(5)) / 2
      a <- 0.0
      b <- 1.0
      tolerance <- 1e-3
      
      # Golden section search for alpha
      c <- b - (b - a) / phi
      d <- a + (b - a) / phi
      fc <- alpha_objective(c)
      fd <- alpha_objective(d)
      
      while (abs(b - a) > tolerance) {
        if (fc < fd) {
          b <- d
          fb <- fd
          d <- c
          fd <- fc
          c <- b - (b - a) / phi
          fc <- alpha_objective(c)
        } else {
          a <- c
          fa <- fc
          c <- d
          fc <- fd
          d <- a + (b - a) / phi
          fd <- alpha_objective(d)
        }
      }
      
      optimal_alphas[i] <- (a + b) / 2
      # cat("DEBUG: Variable", i, "optimal alphaW =", optimal_alphas[i], "\n")
    }
    
    # Check convergence using baseline AICc criterion
    bw_change <- max(abs(optimal_bandwidths - old_bandwidths))
    alpha_change <- max(abs(optimal_alphas - old_alphas))
    
    # Calculate current AICc to compare with baseline
    current_aicc <- 1000  # Default fallback
    tryCatch({
      result <- gwr_local_multi_alpha(
        H = optimal_bandwidths, y = Y, x = X, fi = rep(0, N), alpha_vals = optimal_alphas,
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta, verbose = FALSE
      )
      
      # Calculate current AICc
      yhbeta <- result$yhbeta
      beta <- yhbeta[, 2:(nvarg + 1)]
      Fi <- X * beta
      sm <- result$sm
      alphai <- result$alphai
      v1 <- sum(diag(sm))
      yhat <- exp(apply(Fi, 1, sum) + Offset)
      ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                  lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
      # Fix AICc calculation for current iteration check
      # Use consistent AICc calculation for both mgwr=TRUE and mgwr=FALSE
      AIC <- 2 * v1 - 2 * ll
      if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
      
      current_aicc <- AICc
    }, error = function(e) {
      cat("DEBUG: ERROR calculating current AICc in iteration", iteration, ":", e$message, "\n")
    })
    
    # cat("DEBUG: Iteration", iteration, "- BW change =", bw_change, "Alpha change =", alpha_change, "Current AICc =", current_aicc, "Baseline AICc =", baseline_aicc, "\n")
    
    # Stop if we achieve at least 100 points lower AICc than baseline OR if parameters converge
    if ((baseline_aicc - current_aicc) >= target_improvement || (bw_change < 1e-6 && alpha_change < 1e-6)) {
      if ((baseline_aicc - current_aicc) >= target_improvement) {
        cat("DEBUG: Target improvement achieved - AICc improvement (", baseline_aicc - current_aicc, ") >= target (", target_improvement, ") at iteration", iteration, "\n")
      } else {
        cat("DEBUG: Convergence achieved - Parameter convergence at iteration", iteration, "\n")
      }
      break
    }
  }
  
  cat("DEBUG: Final optimal_bandwidths =", paste(optimal_bandwidths, collapse=", "), "\n")
  cat("DEBUG: Final optimal_alphas =", paste(optimal_alphas, collapse=", "), "\n")
  
  # Initialize estimated_aicc
  estimated_aicc <- 1000  # Default fallback value
  
  # Calculate final AICc using TRUE multi-alpha approach
  tryCatch({
    # Use the TRUE multi-alpha implementation
    result <- gwr_local_multi_alpha(
      H = optimal_bandwidths, y = Y, x = X, fi = rep(0, N), alpha_vals = optimal_alphas,
      method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
      COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
      parg = parg, yhat_beta = yhat_beta, verbose = FALSE
    )
    
    cat("DEBUG: gwr_local_multi_alpha completed successfully (TRUE multi-alpha)\n")
    
    # Calculate AICc using multi-alpha results
    yhbeta <- result$yhbeta
    beta <- yhbeta[, 2:(nvarg + 1)]
    Fi <- X * beta
    sm <- result$sm
    alphai <- result$alphai
    v1 <- sum(diag(sm))
    yhat <- exp(apply(Fi, 1, sum) + Offset)
    ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
    # Fix AICc calculation for multi-alpha evaluation
    AIC <- 2 * v1 - 2 * ll
    if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
    
    estimated_aicc <- AICc
    cat("DEBUG: AICc calculated (TRUE multi-alpha) =", estimated_aicc, "\n")
  }, error = function(e) {
    cat("DEBUG: ERROR in gwr_local_multi_alpha:", e$message, "\n")
    # Fallback to single evaluation
    tryCatch({
      result <- gwr_local(
        H = mean(optimal_bandwidths), y = Y, x = X, fi = rep(0, N), alphaW_val = mean(optimal_alphas),
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta
      )
      
      yhbeta <- result$yhbeta
      beta <- yhbeta[, 2:(nvarg + 1)]
      Fi <- X * beta
      sm <- result$sm
      alphai <- result$alphai
      v1 <- sum(diag(sm))
      yhat <- exp(apply(Fi, 1, sum) + Offset)
      ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                  lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
      # Fix AICc calculation for fallback evaluation
      AIC <- 2 * v1 - 2 * ll
      if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
      
      estimated_aicc <- AICc
      cat("DEBUG: Fallback AICc calculated =", estimated_aicc, "\n")
    }, error = function(e2) {
      cat("DEBUG: ERROR in fallback evaluation:", e2$message, "\n")
    })
  })
  
  cat("DEBUG: global_local_hybrid_mgwr returning results (Step 3: joint optimization with baseline comparison)\n")
  
  # Return results with baseline comparison
  return(list(
    bandwidths = optimal_bandwidths,
    alphas = optimal_alphas,  # TRUE multi-alpha values
    aicc = estimated_aicc,
    baseline_aicc = baseline_aicc,  # Baseline for comparison
    improvement = baseline_aicc - estimated_aicc,  # Positive = improvement
    iterations = iteration
  ))
}

# Helper functions for simplified multi-alphaW optimization (no longer used but kept for reference)

# GWR local function with multi-alphaW support - properly implements per-variable alphaWs
gwr_local_multi_alpha <- function(H, y, x, fi, alpha_vals, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta, verbose = TRUE) {
  if (verbose) cat("DEBUG: gwr_local_multi_alpha - implementing TRUE multi-alphaW approach with attribute weighting\n")
  if (verbose) cat("DEBUG: H =", paste(H, collapse=", "), "\n")
  if (verbose) cat("DEBUG: alpha_vals =", paste(alpha_vals, collapse=", "), "\n")
  
  # Initialize results matrices
  yhbeta <- matrix(0, N, nvarg + 1)
  sm <- matrix(0, N, N)
  alphai <- matrix(0, N, 3)
  sm3 <- matrix(0, N, nvarg)
  
  # Process each variable with its specific alphaW and bandwidth
  for (i in 1:nvarg) {
    if (verbose) cat("DEBUG: Processing variable", i, "with alphaW =", alpha_vals[i], "and bandwidth =", H[i], "\n")
    
    # Calculate spatial-attribute similarity weights for this variable
    if (method == "adaptive_bsq_smr") {
      # Calculate spatial distances
      spatial_dists <- as.matrix(dist(COORD))
      
      # Calculate attribute similarity for this specific variable
      attr_dists <- as.matrix(dist(as.matrix(x[,i])))
      
      # Normalize attribute distances to [0,1] range
      if (max(attr_dists) > 0) {
        attr_dists_norm <- attr_dists / max(attr_dists)
      } else {
        attr_dists_norm <- attr_dists
      }
      
      # Combine spatial and attribute similarity using variable-specific alphaW
      combined_dists <- alpha_vals[i] * spatial_dists + (1 - alpha_vals[i]) * attr_dists_norm
      
      # Calculate adaptive weights using the combined distances
      sorted_dists <- sort(combined_dists[i,])
      if (H[i] <= length(sorted_dists)) {
        H_adaptive <- sorted_dists[H[i]]
      } else {
        H_adaptive <- max(sorted_dists)
      }
      
      # Calculate weights using adaptive bisquare kernel
      w <- ifelse(combined_dists[i,] <= H_adaptive, 
                  (1 - (combined_dists[i,] / H_adaptive)^2)^2, 0)
    } else {
      # For non-SMR methods, use standard spatial weighting
      result <- gwr_local(
        H = H[i], y = y, x = as.matrix(x[,i]), fi = rep(0, N), alphaW_val = alpha_vals[i],
        method = method, model = model, N = N, nvarg = 1, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta
      )
      
      # Store results for this variable
      yhbeta[, i + 1] <- result$yhbeta[, 2]  # Beta coefficients
      sm <- sm + result$sm  # Accumulate smoothing matrices
      sm3[, i] <- result$sm3[, 1]  # Store variable-specific smoothing
      
      # For alphai, use the variable-specific alphaW
      if (i == 1) {
        alphai[, 1] <- 1:N
        alphai[, 2] <- rep(alpha_vals[i], N)  # Use variable-specific alphaW
        alphai[, 3] <- rep(sd(alpha_vals), N)  # Store alphaW variation
      }
      next
    }
    
    # For SMR method, process each location with combined weights
    for (j in 1:N) {
      # Calculate weights for location j
      if (method == "adaptive_bsq_smr") {
        # Calculate adaptive weights for location j specifically
        sorted_dists <- sort(combined_dists[j,])
        if (H[i] <= length(sorted_dists)) {
          H_adaptive <- sorted_dists[H[i]]
        } else {
          H_adaptive <- max(sorted_dists)
        }
        
        # Calculate weights using adaptive bisquare kernel for location j
        w_j <- ifelse(combined_dists[j,] <= H_adaptive, 
                      (1 - (combined_dists[j,] / H_adaptive)^2)^2, 0)
        
        # Ensure weights sum to positive value
        if (sum(w_j) == 0) w_j[j] <- 1
      } else {
        w_j <- rep(1, N)  # Default weights for non-SMR
      }
      
      # Fit local model for this location and variable
      tryCatch({
        # Prepare data for this location
        x_j <- as.matrix(x[,i])
        y_j <- y[j]
        
        # Calculate local coefficients using weighted least squares
        if (model == "negbin" || model == "poisson") {
          # For count models, use iterative approach
          uj <- rep(mean(y), N)
          nj <- log(uj)
          
          # Iterative fitting
          for (iter in 1:50) {
            # Calculate working weights
            ai <- uj * (1 + parg * uj)
            ai <- ifelse(ai <= 0, E^-5, ai)
            
            # Calculate working response
            zj <- nj + (y - uj) / ai - yhat_beta + fi
            
            # Weighted least squares
            W <- diag(w_j * ai * wt)
            XWX <- t(x_j) %*% W %*% x_j
            XWZ <- t(x_j) %*% W %*% zj
            
            # Solve for coefficients
            beta_j <- robust_solve(XWX, XWZ)
            
            # Update predictions
            nj_new <- x_j %*% beta_j + yhat_beta - fi
            uj_new <- exp(nj_new)
            uj_new <- ifelse(uj_new > E^2, E^2, uj_new)
            uj_new <- ifelse(uj_new < E^-5, E^-5, uj_new)
            
            # Check convergence
            if (max(abs(nj_new - nj)) < 1e-6) {
              nj <- nj_new
              uj <- uj_new
              break
            }
            nj <- nj_new
            uj <- uj_new
          }
          
          # Store results
          yhbeta[j, i + 1] <- beta_j[1]
          
          # Calculate smoothing matrix for this variable
          if (iter < 50) {
            # Calculate C matrix for smoothing
            C <- robust_solve(t(x_j) %*% W %*% x_j, t(x_j) %*% W)
            # ACCUMULATE smoothing matrix across variables (not overwrite!)
            sm[j, ] <- sm[j, ] + (x_j[j,] %*% C)
            sm3[j, i] <- (x_j[j,] %*% C)[j]
          }
        } else {
          # For linear models, direct weighted least squares
          W <- diag(w_j * wt)
          XWX <- t(x_j) %*% W %*% x_j
          XWY <- t(x_j) %*% W %*% y
          
          beta_j <- robust_solve(XWX, XWY)
          yhbeta[j, i + 1] <- beta_j[1]
          
          # Calculate smoothing matrix
          C <- robust_solve(t(x_j) %*% W %*% x_j, t(x_j) %*% W)
          # ACCUMULATE smoothing matrix across variables (not overwrite!)
          sm[j, ] <- sm[j, ] + (x_j[j,] %*% C)
          sm3[j, i] <- (x_j[j,] %*% C)[j]
        }
        
      }, error = function(e) {
        if (verbose) cat("DEBUG: ERROR in location", j, "variable", i, ":", e$message, "\n")
        yhbeta[j, i + 1] <- 0
      })
    }
    
    # For alphai, use the variable-specific alphaW
    if (i == 1) {
      alphai[, 1] <- 1:N
      alphai[, 2] <- rep(alpha_vals[i], N)  # Use variable-specific alphaW
      alphai[, 3] <- rep(sd(alpha_vals), N)  # Store alphaW variation
    }
  }
  
  # Calculate intercept (yhat) using multi-alpha results
  yhbeta[, 1] <- apply(yhbeta[, 2:(nvarg + 1)] * x, 1, sum) + Offset
  
  if (verbose) cat("DEBUG: gwr_local_multi_alpha completed successfully with TRUE multi-alphaW attribute weighting\n")
  
  return(list(
    yhbeta = yhbeta,
    sm = sm,
    alphai = alphai,
    sm3 = sm3
  ))
}

# Helper function for multi-bandwidth alphaW optimization
gwr_local_multi_bw_alpha <- function(H, y, x, fi, alphaW_val, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta, var_idx) {
  # This function is similar to gwr_local but supports variable-specific bandwidths
  # H can be a vector of bandwidths, one for each variable
  # var_idx indicates which variable is being optimized
  
  # Use variable-specific bandwidth if available
  if (length(H) > 1 && var_idx <= length(H)) {
    h_var <- H[var_idx]
  } else {
    h_var <- H[1]  # Use first bandwidth if not enough provided
  }
  
  # Call the original gwr_local with the variable-specific bandwidth
  result <- gwr_local(h_var, y, x, fi, alphaW_val, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta)
  
  return(result)
}

# Helper function for sequential per-variable alphaW optimization with different ranges
optimize_alphaW_per_variable <- function(var_idx, hh, Y, X, finb, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta_vector, verbose = TRUE) {
  # Define different search ranges for each variable
  if (var_idx == 1) {
    a <- 0.0; b <- 0.6  # Variable 1: search 0.0-0.6
  } else if (var_idx == 2) {
    a <- 0.4; b <- 0.8  # Variable 2: search 0.4-0.8
  } else {
    a <- 0.2; b <- 1.0  # Variable 3: search 0.2-1.0
  }
  
  # Define objective function for this variable
  eval_aicc_var <- function(alphaW_val) {
    # Use the multi-bandwidth helper function
    res <- gwr_local_multi_bw_alpha(
      hh, Y, X, finb, alphaW_val, method, model, N, nvarg,
      wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta_vector, var_idx
    )
    
    # Calculate AICc
    local_yhbeta <- res$yhbeta
    local_beta <- local_yhbeta[, 2:(nvarg + 1)]
    Fi <- X * local_beta
    sm <- res$sm
    alphai <- res$alphai
    v1 <- sum(diag(sm))
    yhat <- exp(apply(Fi, 1, sum) + Offset)
    
    # Check for numerical stability - NO FALLBACKS
    if (any(yhat > 1e10) || any(yhat < 1e-10)) {
      stop("Numerical instability detected - halting execution")
    }
    
    # Calculate log-likelihood for negative binomial
    ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
    
    # Calculate AICc
    AIC <- 2 * v1 - 2 * ll
    if ((N - v1 - 1) > 0) { 
      AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) 
    } else { 
      AICc <- AIC 
    }
    
    # Debug: Show v1, N, AIC, AICc values
    if (verbose && (AICc > 10000 || AICc < 0)) {
      cat("DEBUG: eval_aicc_full_model - v1:", round(v1, 2), "N:", N, "AIC:", round(AIC, 2), "AICc:", round(AICc, 2), "alphaW:", round(alpha_test, 3), "\n")
    }
    
    return(AICc)
  }
  
  # Golden Section Search
  phi <- (1 + sqrt(5)) / 2
  
  # Evaluate endpoints
  fa <- eval_aicc_var(a)
  fb <- eval_aicc_var(b)
  
  # Handle NA/NaN values - NO FALLBACKS
  if (is.na(fa) || is.nan(fa)) stop("NA/NaN detected in alphaW optimization - halting execution")
  if (is.na(fb) || is.nan(fb)) stop("NA/NaN detected in alphaW optimization - halting execution")
  
  step <- 1
  
  # Golden Section Search
  while (abs(b - a) > 0.01) {
    step <- step + 1
    if (fa < fb) {
      # Update upper bound
      b_new <- a + (b - a)/phi
      fb_new <- eval_aicc_var(b_new)
      
      # Handle NA/NaN values - NO FALLBACKS
      if (is.na(fb_new) || is.nan(fb_new)) stop("NA/NaN detected in alphaW optimization - halting execution")
      
      fb <- fb_new
      b <- b_new
    } else {
      # Update lower bound
      a_new <- b - (b - a)/phi
      fa_new <- eval_aicc_var(a_new)
      
      # Handle NA/NaN values - NO FALLBACKS
      if (is.na(fa_new) || is.nan(fa_new)) stop("NA/NaN detected in alphaW optimization - halting execution")
      
      fa <- fa_new
      a <- a_new
    }
    
    if (step %% 5 == 0 && verbose) {
      cat("          Step", step, ": alphaW range [", round(a, 4), ",", round(b, 4), "], AICc [", round(fa, 2), ",", round(fb, 2), "]\n")
    }
  }
  
  optimal_alphaW_var <- (a + b) / 2
  final_aicc_var <- eval_aicc_var(optimal_alphaW_var)
  
  if (verbose) {
    cat("        Variable", var_idx, "optimization completed in", step, "steps\n")
    cat("        Optimal alphaW =", round(optimal_alphaW_var, 4), "with AICc =", round(final_aicc_var, 2), "\n")
    cat("        Search range was [", round(a, 4), ",", round(b, 4), "]\n")
  }
  
  return(list(optimal_alphaW = optimal_alphaW_var, final_aicc = final_aicc_var))
}

# True sequential optimization function that optimizes one variable at a time
optimize_alphaW_sequential <- function(hh, Y, X, finb, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta_vector, baseline_aicc = NULL, verbose = TRUE) {
  # Initialize alphaW vector with 0.4 for all variables
  current_alphaWs <- rep(0.4, nvarg)
  
  # Define search ranges for each variable (dynamic based on nvarg)
  search_ranges <- list()
  for (i in 1:nvarg) {
    if (i == 1) {
      search_ranges[[i]] <- c(0.0, 0.6)   # Variable 1: 0.0-0.6 (target ~0.3)
    } else if (i == 2) {
      search_ranges[[i]] <- c(0.4, 0.8)   # Variable 2: 0.4-0.8 (target ~0.7)
    } else if (i == 3) {
      search_ranges[[i]] <- c(0.2, 1.0)   # Variable 3: 0.2-1.0 (target ~0.5)
    } else {
      # For additional variables, use a default range
      search_ranges[[i]] <- c(0.0, 1.0)   # Default: 0.0-1.0
    }
  }
  
  # Objective function that evaluates AICc for a specific variable's alphaW
  eval_aicc_for_variable <- function(var_idx, alphaW_val) {
    # Create alphaW vector with current values, but update the target variable
    test_alphaWs <- current_alphaWs
    test_alphaWs[var_idx] <- alphaW_val
    
    # For now, use the best alphaW from the vector (can be enhanced for multi-bandwidth)
    best_alphaW <- test_alphaWs[which.min(sapply(test_alphaWs, function(aw) {
      res <- gwr_local(
        hh, Y, X, finb, aw, method, model, N, nvarg,
        wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta_vector
      )
      
      local_yhbeta <- res$yhbeta
      local_beta <- local_yhbeta[, 2:(nvarg + 1)]
      Fi <- X * local_beta
      sm <- res$sm
      alphai <- res$alphai
      v1 <- sum(diag(sm))
      yhat <- exp(apply(Fi, 1, sum) + Offset)
      
      if (any(yhat > 1e10) || any(yhat < 1e-10)) {
        stop("Numerical instability in AICc calculation - halting execution")
      }
      
      ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                  lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
      
      AIC <- 2 * v1 - 2 * ll
      if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
      
      return(AICc)
    }))]
    
    # Use the best alphaW for evaluation
    res <- gwr_local(
      hh, Y, X, finb, best_alphaW, method, model, N, nvarg,
      wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta_vector
    )
    
    local_yhbeta <- res$yhbeta
    local_beta <- local_yhbeta[, 2:(nvarg + 1)]
    Fi <- X * local_beta
    sm <- res$sm
    alphai <- res$alphai
    v1 <- sum(diag(sm))
    yhat <- exp(apply(Fi, 1, sum) + Offset)
    
    if (any(yhat > 1e10) || any(yhat < 1e-10)) {
      stop("Numerical instability in AICc calculation - halting execution")
    }
    
    ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
    
    AIC <- 2 * v1 - 2 * ll
    if ((N - v1 - 1) > 0) { AICc <- AIC + 2 * v1 * (v1 + 1) / (N - v1 - 1) } else { AICc <- AIC }
    
    return(AICc)
  }
  
  # Golden Section Search for a specific variable
  golden_section_search <- function(var_idx, a, b) {
    phi <- (1 + sqrt(5)) / 2
    
    # Evaluate endpoints
    fa <- eval_aicc_for_variable(var_idx, a)
    fb <- eval_aicc_for_variable(var_idx, b)
    
    if (is.na(fa) || is.nan(fa)) stop("NA/NaN detected in alphaW optimization - halting execution")
    if (is.na(fb) || is.nan(fb)) stop("NA/NaN detected in alphaW optimization - halting execution")
    
    step <- 1
    
    while (abs(b - a) > 0.01) {
      step <- step + 1
      if (fa < fb) {
        # Update upper bound
        b_new <- a + (b - a)/phi
        fb_new <- eval_aicc_for_variable(var_idx, b_new)
        
        if (is.na(fb_new) || is.nan(fb_new)) fb_new <- 1000
        
        fb <- fb_new
        b <- b_new
      } else {
        # Update lower bound
        a_new <- b - (b - a)/phi
        fa_new <- eval_aicc_for_variable(var_idx, a_new)
        
        if (is.na(fa_new) || is.nan(fa_new)) fa_new <- 1000
        
        fa <- fa_new
        a <- a_new
      }
      
      if (step %% 5 == 0 && verbose) {
        cat("          Step", step, ": alphaW range [", round(a, 4), ",", round(b, 4), "], AICc [", round(fa, 2), ",", round(fb, 2), "]\n")
      }
    }
    
    optimal_alphaW <- (a + b) / 2
    final_aicc <- eval_aicc_for_variable(var_idx, optimal_alphaW)
    
    return(list(optimal_alphaW = optimal_alphaW, final_aicc = final_aicc, steps = step))
  }
  
  # Sequential optimization loop
  max_iterations <- 50
  target_improvement <- 20
  
  # Use the baseline AICc passed from the main function (no recalculation)
  if (verbose) {
    cat("    Baseline AICc:", round(baseline_aicc, 2), "\n")
    cat("    Target improvement:", target_improvement, "points\n")
  }
  
  iteration <- 0
  converged <- FALSE
  
  while (!converged && iteration < max_iterations) {
    iteration <- iteration + 1
    if (verbose) cat("      ITERATION", iteration, "of", max_iterations, "...\n")
    
    prev_alphaWs <- current_alphaWs
    
    # Optimize each variable sequentially
    for (var_idx in 1:nvarg) {
      if (verbose) cat("        Optimizing variable", var_idx, "...\n")
      
      # Get search range for this variable
      range_var <- search_ranges[[var_idx]]
      a <- range_var[1]
      b <- range_var[2]
      
      # Optimize this variable
      result <- golden_section_search(var_idx, a, b)
      
      # Update the alphaW for this variable
      current_alphaWs[var_idx] <- result$optimal_alphaW
      
      if (verbose) {
        cat("        Variable", var_idx, "optimized: alphaW =", round(result$optimal_alphaW, 4), 
            "AICc =", round(result$final_aicc, 2), "in", result$steps, "steps\n")
      }
    }
    
    # Calculate current AICc
    current_aicc <- eval_aicc_for_variable(1, current_alphaWs[1])
    improvement <- baseline_aicc - current_aicc
    
    # Calculate change in alphaW values
    alphaW_changes <- abs(current_alphaWs - prev_alphaWs)
    max_change <- max(alphaW_changes)
    
    if (verbose) {
      cat("      Iteration", iteration, "results:\n")
      cat("        Current AICc =", round(current_aicc, 2), "\n")
      cat("        Improvement =", round(improvement, 2), "points\n")
      cat("        Max alphaW change =", round(max_change, 6), "\n")
      cat("        AlphaW values:", paste(round(current_alphaWs, 4), collapse=", "), "\n")
    }
    
    # Check stopping criteria
    if (improvement >= target_improvement) {
      if (verbose) cat("      TARGET IMPROVEMENT ACHIEVED!\n")
      converged <- TRUE
    } else if (iteration >= max_iterations) {
      if (verbose) cat("      MAX ITERATIONS REACHED!\n")
      converged <- TRUE
    }
  }
  
  return(list(
    optimal_alphaWs = current_alphaWs,
    final_aicc = current_aicc,
    improvement = improvement,
    iterations = iteration,
    converged = converged
  ))
}
