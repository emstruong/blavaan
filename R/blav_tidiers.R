# Tidy and glance methods for blavaan objects
# Following broom ecosystem conventions for Bayesian models

# Define tidy generic if not already available (e.g., from broom or generics)
# This ensures the S3 method registration works even without those packages
#' @export
tidy <- function(x, ...) UseMethod("tidy")

#' @export
glance <- function(x, ...) UseMethod("glance")

#' Tidy a blavaan object
#'
#' @param x A \code{blavaan} object
#' @param conf.int Logical indicating whether to include credible intervals.
#'   Default is \code{TRUE}.
#' @param conf.level The probability mass to include in the credible interval.
#'   Default is 0.95.
#' @param conf.method Method for computing credible intervals. Either "quantile"
#'   (equal-tailed interval) or "HPDinterval" (highest posterior density).
#'   Default is "quantile".
#' @param standardized Logical indicating whether to include standardized
#'   estimates. Default is \code{FALSE}.
#' @param rhat Logical indicating whether to include the Rhat convergence
#'   diagnostic. Default is \code{TRUE}.
#' @param ess Logical indicating whether to include the effective sample size.
#'   Default is \code{TRUE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{tibble} (if available) or \code{data.frame} with columns:
#'   \item{term}{Parameter name (lhs, op, rhs combined)}
#'   \item{op}{Operator from the model syntax}
#'   \item{estimate}{Posterior mean}
#'   \item{std.error}{Posterior standard deviation}
#'   \item{conf.low}{Lower bound of credible interval (if \code{conf.int = TRUE})}
#'   \item{conf.high}{Upper bound of credible interval (if \code{conf.int = TRUE})}
#'   \item{std.all}{Standardized estimate (if \code{standardized = TRUE})}
#'   \item{rhat}{Rhat convergence diagnostic (if \code{rhat = TRUE})}
#'   \item{ess}{Effective sample size (if \code{ess = TRUE})}
#'   \item{prior}{Prior distribution specification}
#'   \item{group}{Group number (for multigroup models)}
#'
#' @examples
#' \dontrun{
#' library(blavaan)
#'
#' HS.model <- 'visual =~ x1 + x2 + x3'
#' fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123)
#' tidy(fit)
#' tidy(fit, conf.int = TRUE, conf.level = 0.90)
#' tidy(fit, standardized = TRUE)
#' }
#'
#' @export
tidy.blavaan <- function(x, conf.int = TRUE, conf.level = 0.95,
                         conf.method = c("quantile", "HPDinterval"),
                         standardized = FALSE, rhat = TRUE, ess = TRUE, ...) {

  conf.method <- match.arg(conf.method)

  # Get parameter estimates
  PE <- parameterEstimates(x, se = TRUE, zstat = FALSE, pvalue = FALSE,
                           ci = FALSE, standardized = standardized,
                           remove.eq = FALSE, remove.system.eq = TRUE,
                           remove.ineq = FALSE, remove.def = FALSE)

  # Build the base data frame
  result <- data.frame(
    term = paste0(PE$lhs, PE$op, PE$rhs),
    op = PE$op,
    estimate = PE$est,
    std.error = PE$se,
    stringsAsFactors = FALSE
  )

  # Add credible intervals if requested
  if (conf.int) {
    ci_bounds <- compute_blavaan_ci(x, conf.level, conf.method)
    if (!is.null(ci_bounds)) {
      result$conf.low <- ci_bounds$lower
      result$conf.high <- ci_bounds$upper
    } else {
      result$conf.low <- NA_real_
      result$conf.high <- NA_real_
    }
  }

  # Add standardized estimates if requested
  if (standardized && "std.all" %in% names(PE)) {
    result$std.all <- PE$std.all
  }

  # Add Rhat if requested
  if (rhat) {
    result$rhat <- get_blavaan_rhat(x, PE)
  }

  # Add effective sample size if requested
  if (ess) {
    result$ess <- get_blavaan_ess(x, PE)
  }

  # Add prior information
  result$prior <- get_blavaan_priors(x, PE)

  # Add group information for multigroup models
  if ("group" %in% names(PE)) {
    result$group <- PE$group
  }

  # Convert to tibble if available
  as_blavaan_tibble(result)
}


#' Glance at a blavaan object
#'
#' @param x A \code{blavaan} object
#' @param fit.indices Character vector of fit indices to compute from
#'   \code{blavFitIndices}. Use \code{"none"} to skip Bayesian fit indices
#'   (faster). Use \code{"default"} for BRMSEA, BGammaHat, and BMc.
#'   Default is \code{"none"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A single-row \code{tibble} (if available) or \code{data.frame}
#'   with columns:
#'   \item{npar}{Number of estimated parameters}
#'   \item{nobs}{Total number of observations}
#'   \item{ngroups}{Number of groups}
#'   \item{estimator}{Estimation method used}
#'   \item{ppp}{Posterior predictive p-value}
#'   \item{looic}{Leave-one-out information criterion}
#'   \item{p_loo}{Effective number of parameters (LOO)}
#'   \item{waic}{Widely applicable information criterion}
#'   \item{p_waic}{Effective number of parameters (WAIC)}
#'   \item{dic}{Deviance information criterion}
#'   \item{p_dic}{Effective number of parameters (DIC)}
#'   \item{bic}{Bayesian information criterion}
#'   \item{margloglik}{Marginal log-likelihood}
#'   \item{converged}{Logical indicating convergence (all Rhat < 1.05)}
#'   \item{nchains}{Number of MCMC chains}
#'   Additional columns from \code{blavFitIndices} are included if requested.
#'
#' @examples
#' \dontrun{
#' library(blavaan)
#'
#' HS.model <- 'visual =~ x1 + x2 + x3'
#' fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123)
#' glance(fit)
#' glance(fit, fit.indices = "default")  # includes BRMSEA etc.
#' }
#'
#' @export
glance.blavaan <- function(x, fit.indices = "none", ...) {

  # Get basic model information
  ngroups <- blavInspect(x, "ngroups")
  nobs <- blavInspect(x, "ntotal")
  nchains <- blavInspect(x, "n.chains")

  # Get fit measures
  bopts <- blavInspect(x, "options")
  test_available <- bopts$test != "none"

  # Initialize result
  result <- data.frame(
    npar = as.integer(lavInspect(x, "npar")),
    nobs = as.integer(nobs),
    ngroups = as.integer(ngroups),
    estimator = "Bayes",
    stringsAsFactors = FALSE
  )

  # Get available fit measures
  if (test_available) {
    fm <- tryCatch(
      fitMeasures(x, fit.measures = "all"),
      error = function(e) NULL
    )

    if (!is.null(fm)) {
      # Posterior predictive p-value
      if ("ppp" %in% names(fm)) {
        result$ppp <- unname(fm["ppp"])
      }

      # Information criteria
      if ("looic" %in% names(fm)) {
        result$looic <- unname(fm["looic"])
        result$p_loo <- unname(fm["p_loo"])
      }
      if ("waic" %in% names(fm)) {
        result$waic <- unname(fm["waic"])
        result$p_waic <- unname(fm["p_waic"])
      }
      if ("dic" %in% names(fm)) {
        result$dic <- unname(fm["dic"])
        result$p_dic <- unname(fm["p_dic"])
      }
      if ("bic" %in% names(fm)) {
        result$bic <- unname(fm["bic"])
      }
      if ("margloglik" %in% names(fm)) {
        result$margloglik <- unname(fm["margloglik"])
      }
    }
  }

  # Check convergence (Rhat < 1.05 for all parameters)
  rhat_vals <- tryCatch(
    blavInspect(x, "rhat"),
    error = function(e) NULL
  )
  if (!is.null(rhat_vals)) {
    result$converged <- all(rhat_vals < 1.05, na.rm = TRUE)
  }

  result$nchains <- as.integer(nchains)

  # Add Bayesian fit indices if requested
  if (!identical(fit.indices, "none")) {
    bfi <- tryCatch({
      if (identical(fit.indices, "default")) {
        blavFitIndices(x, fit.measures = c("BRMSEA", "BGammaHat", "BMc"))
      } else {
        blavFitIndices(x, fit.measures = fit.indices)
      }
    }, error = function(e) NULL)

    if (!is.null(bfi)) {
      bfi_summary <- summary(bfi, central.tendency = "mean", hpd = FALSE)
      # Add EAP (posterior mean) for each fit index
      for (i in seq_len(nrow(bfi_summary))) {
        idx_name <- rownames(bfi_summary)[i]
        result[[idx_name]] <- bfi_summary[i, "EAP"]
      }
    }
  }

  as_blavaan_tibble(result)
}


# Helper function to compute credible intervals
compute_blavaan_ci <- function(x, conf.level = 0.95, method = "quantile") {

  jagtarget <- lavInspect(x, "options")$target == "jags"
  newpt <- x@ParTable
  if (!("group" %in% names(newpt))) newpt$group <- rep(1, length(newpt$lhs))
  if (!("level" %in% names(newpt))) newpt$level <- rep("within", length(newpt$lhs))
  newpt$group[newpt$group == 0] <- 1

  # Get parameter entries
  if (jagtarget) {
    pte2 <- which(!is.na(newpt$jagpnum))
  } else {
    pte2 <- which(newpt$free > 0)
  }

  # Get parameter estimates for matching
  PE <- parameterEstimates(x, se = FALSE, zstat = FALSE, pvalue = FALSE,
                           ci = FALSE, remove.eq = FALSE, remove.system.eq = TRUE,
                           remove.ineq = FALSE, remove.def = FALSE)
  if (!("group" %in% names(PE))) PE$group <- 1
  if (!("level" %in% names(PE))) PE$level <- "within"
  PE$group[PE$group == 0] <- 1

  peentry <- match(
    paste(newpt$lhs[pte2], newpt$op[pte2], newpt$rhs[pte2],
          newpt$group[pte2], newpt$level[pte2], sep = ""),
    paste(PE$lhs, PE$op, PE$rhs, PE$group, PE$level, sep = "")
  )

  # Initialize CI vectors
  ci_lower <- rep(NA_real_, nrow(PE))
  ci_upper <- rep(NA_real_, nrow(PE))

  if (method == "HPDinterval") {
    # Use HPD intervals
    hpd <- tryCatch(
      blavInspect(x, "hpd", level = conf.level),
      error = function(e) NULL
    )
    if (!is.null(hpd)) {
      ci_lower[peentry] <- hpd[, "lower"]
      ci_upper[peentry] <- hpd[, "upper"]
    }
  } else {
    # Use quantile intervals
    alpha <- (1 - conf.level) / 2
    probs <- c(alpha, 1 - alpha)

    if (jagtarget) {
      if (!is.null(x@external$mcmcout$HPD)) {
        # JAGS stores 95% intervals
        if (conf.level == 0.95 && 'Lower95' %in% colnames(x@external$mcmcout$HPD)) {
          ci_lower[peentry] <- x@external$mcmcout$HPD[newpt$jagpnum[pte2], 'Lower95']
          ci_upper[peentry] <- x@external$mcmcout$HPD[newpt$jagpnum[pte2], 'Upper95']
        } else {
          # Compute from MCMC samples
          draws <- blavInspect(x, "mcmc")
          draws <- do.call("rbind", draws)
          ci_bounds <- t(apply(draws, 2, quantile, probs = probs))
          ci_lower[peentry] <- ci_bounds[, 1]
          ci_upper[peentry] <- ci_bounds[, 2]
        }
      }
    } else {
      # Stan
      parsumm <- rstan::summary(x@external$mcmcout)
      pct_cols <- paste0(probs * 100, "%")
      if (all(pct_cols %in% colnames(parsumm$summary))) {
        ci_lower[peentry] <- parsumm$summary[newpt$stansumnum[pte2], pct_cols[1]]
        ci_upper[peentry] <- parsumm$summary[newpt$stansumnum[pte2], pct_cols[2]]
      } else if (conf.level == 0.95 && all(c("2.5%", "97.5%") %in% colnames(parsumm$summary))) {
        ci_lower[peentry] <- parsumm$summary[newpt$stansumnum[pte2], "2.5%"]
        ci_upper[peentry] <- parsumm$summary[newpt$stansumnum[pte2], "97.5%"]
      } else {
        # Compute from MCMC samples
        draws <- blavInspect(x, "mcmc")
        draws <- do.call("rbind", draws)
        ci_bounds <- t(apply(draws, 2, quantile, probs = probs))
        ci_lower[peentry] <- ci_bounds[, 1]
        ci_upper[peentry] <- ci_bounds[, 2]
      }
    }
  }

  list(lower = ci_lower, upper = ci_upper)
}


# Helper function to get Rhat values matched to parameter estimates
get_blavaan_rhat <- function(x, PE) {
  rhat_vals <- tryCatch(
    blavInspect(x, "rhat"),
    error = function(e) NULL
  )

  if (is.null(rhat_vals)) {
    return(rep(NA_real_, nrow(PE)))
  }

  rhat_result <- rep(NA_real_, nrow(PE))
  rhat_names <- names(rhat_vals)

  # Match by parameter label
  pe_labels <- paste0(PE$lhs, PE$op, PE$rhs)
  matched_idx <- match(pe_labels, rhat_names)
  rhat_result[!is.na(matched_idx)] <- rhat_vals[matched_idx[!is.na(matched_idx)]]

  rhat_result
}


# Helper function to get effective sample size matched to parameter estimates
get_blavaan_ess <- function(x, PE) {
  ess_vals <- tryCatch(
    blavInspect(x, "neff"),
    error = function(e) NULL
  )

  if (is.null(ess_vals)) {
    return(rep(NA_real_, nrow(PE)))
  }

  ess_result <- rep(NA_real_, nrow(PE))
  ess_names <- names(ess_vals)

  # Match by parameter label
  pe_labels <- paste0(PE$lhs, PE$op, PE$rhs)
  matched_idx <- match(pe_labels, ess_names)
  ess_result[!is.na(matched_idx)] <- ess_vals[matched_idx[!is.na(matched_idx)]]

  ess_result
}


# Helper function to get prior specifications matched to parameter estimates
get_blavaan_priors <- function(x, PE) {
  pt <- x@ParTable
  prior_result <- rep(NA_character_, nrow(PE))

  if (!("prior" %in% names(pt))) {
    return(prior_result)
  }

  # Match parameters
  pe_labels <- paste0(PE$lhs, PE$op, PE$rhs)
  pt_labels <- paste0(pt$lhs, pt$op, pt$rhs)

  if ("group" %in% names(PE) && "group" %in% names(pt)) {
    pe_labels <- paste0(pe_labels, ".g", PE$group)
    pt_labels <- paste0(pt_labels, ".g", pt$group)
  }

  matched_idx <- match(pe_labels, pt_labels)
  prior_result[!is.na(matched_idx)] <- pt$prior[matched_idx[!is.na(matched_idx)]]

  # Replace empty strings with NA
  prior_result[prior_result == ""] <- NA_character_

  prior_result
}


# Helper function to convert to tibble if available
as_blavaan_tibble <- function(x) {
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::as_tibble(x)
  } else {
    class(x) <- c("tbl_df", "tbl", "data.frame")
    x
  }
}
