#' ITS analysis for count outcomes with multiple interruptions (no seasonality adjustment)
#'
#' \code{its_poisson_wo_seas} fits a Poisson regression model to an ITS with one or more disruption events. 
#' It returns the fitted model, the summary (including the relative risk), and the original data 
#' appended with the model's factual and counterfactual predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offset term.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL. Default value is NULL.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent (lasting until the end of the series). Default is NULL.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular" (homogenous effects), "multiple" (distinct effects), or "mixed" (hierarchical effects). Default is "singular".
#' @param over_dispersion Logical - indicating whether a quasi-Poisson model should be used to account for overdispersion (TRUE), or not (FALSE). Default value is FALSE. Note: This argument is ignored if multi_model_type is "mixed" due to lme4 limitations.
#' @param impact_model A string specifying the assumed impact model. Options include "full" (level and slope change), "level" (just level change), and "slope" (just slope change). Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed, or just the relevant effect size. Default value is FALSE.
#' @return The function returns a list with three elements: the fitted Poisson regression model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' # Simulating two disruptions
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' lengths <- c(6, 12) 
#' fit <- its_poisson_wo_seas(data=data, form=form, offset_name="labour", time_name="time", 
#'                            intervention_start_indices=start_indices, intervention_lengths=lengths, 
#'                            multi_model_type="multiple", over_dispersion=TRUE, 
#'                            impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export

its_poisson_wo_seas <- function(data, form, offset_name=NULL, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", over_dispersion=FALSE, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  pred <- predC <- indicator <- NULL
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(rlang::is_formula(form))) stop("Please make sure that form is a formula object")
  if (!(time_name %in% colnames(data))) stop("Please make sure that time_name belongs to colnames(data)")
  
  if (!all(data[[time_name]] == (1:nrow(data)))) {
    data[[time_name]] = 1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!is.logical(over_dispersion)) stop("over_dispersion must be either TRUE or FALSE")
  if (!(impact_model %in% c("full","level","slope"))) stop("Please enter a valid impact_model name: 'full', 'level', or 'slope'.")
  
  n <- nrow(data)
  K <- length(intervention_start_indices)
  
  # 1. DATA PREP & FORMULA BUILDER 
  if (multi_model_type %in% c("singular", "multiple")) {
    data <- build_disruption_features(data, time_name, intervention_start_indices, intervention_lengths)
  } else if (multi_model_type == "mixed") {
    data <- build_mixed_features(data, time_name, intervention_start_indices, intervention_lengths)
  }
  
  if (multi_model_type == "singular") {
    data$agg_indicator <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE])
    data$agg_shifted_time <- 0
    for (k in 1:K) data$agg_shifted_time <- data$agg_shifted_time + (data[[paste0("indicator_", k)]] * data[[paste0("shifted_time_", k)]])
    
    if (impact_model == "full") form_update <- update(form, as.formula("~ . + agg_indicator + agg_shifted_time"))
    else if (impact_model == "level") form_update <- update(form, as.formula("~ . + agg_indicator"))
    else if (impact_model == "slope") form_update <- update(form, as.formula("~ . + agg_shifted_time"))
    
  } else if (multi_model_type == "multiple") {
    terms <- c()
    for (k in 1:K) {
      ind <- paste0("indicator_", k)
      shift <- paste0("shifted_time_", k)
      if (impact_model == "full") terms <- c(terms, ind, paste0(ind, ":", shift))
      else if (impact_model == "level") terms <- c(terms, ind)
      else if (impact_model == "slope") terms <- c(terms, paste0(ind, ":", shift))
    }
    form_update <- update(form, as.formula(paste("~ . +", paste(terms, collapse = " + "))))
    
  } else if (multi_model_type == "mixed") {
    if (impact_model == "full") form_update <- update(form, as.formula("~ . + generic_indicator + generic_shifted_time + (0 + generic_indicator + generic_shifted_time | disruption_id)"))
    else if (impact_model == "level") form_update <- update(form, as.formula("~ . + generic_indicator + (0 + generic_indicator | disruption_id)"))
    else if (impact_model == "slope") form_update <- update(form, as.formula("~ . + generic_shifted_time + (0 + generic_shifted_time | disruption_id)"))
  } else stop("Invalid multi_model_type: 'singular', 'multiple', or 'mixed'.")
  
  if (!is.null(offset_name)) form_update <- update(form_update, as.formula(paste0("~ . + offset(log(", sym(offset_name), "))")))
  
  # 2. MODEL FITTING 
  set.seed(1)
  if (multi_model_type == "mixed") {
    if (!requireNamespace("lme4", quietly = TRUE)) stop("The 'lme4' package is required for the mixed model.")
    model <- lme4::glmer(form_update, data, family = poisson)
  } else {
    model <- glm(form_update, data, family = ifelse(!isTRUE(over_dispersion), poisson, quasipoisson))
  }
  
  # 3. EFFECT SIZE EXTRACTION 
  if (multi_model_type %in% c("singular", "multiple")) {
    if (multi_model_type == "singular") data$is_active <- data$agg_indicator > 0
    else data$is_active <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE]) > 0
  } else {
    data$is_active <- data$generic_indicator > 0
  }
  
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  total_L <- sum(intervention_lengths)
  
  if (multi_model_type == "singular") {
    cov_mat <- vcov(model)
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    
    if (impact_model == "full") {
      var_lin <- cov_mat["agg_indicator", "agg_indicator"] + (b_mult^2) * cov_mat["agg_shifted_time", "agg_shifted_time"] + (b_mult * 2) * cov_mat["agg_indicator", "agg_shifted_time"]
      lin_est <- model$coefficients["agg_indicator"] + b_mult * model$coefficients["agg_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat["agg_indicator", "agg_indicator"]
      lin_est <- model$coefficients["agg_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat["agg_shifted_time", "agg_shifted_time"]
      lin_est <- b_mult * model$coefficients["agg_shifted_time"]
    }
    
    RR <- exp(lin_est); CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin)); TS <- lin_est / sqrt(var_lin)
    p_value <- round(2 * pnorm(abs(TS), lower.tail = FALSE), 2)
    ret_vec <- c(RR, exp(CI_lin), p_value); names(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
  } else if (multi_model_type == "multiple") {
    cov_mat <- vcov(model)
    betas <- model$coefficients
    ret_vec <- matrix(NA, nrow = K, ncol = 4)
    rownames(ret_vec) <- paste0("Disruption_", 1:K)
    colnames(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      ind_name <- paste0("indicator_", k)
      shift_name <- paste0("shifted_time_", k)
      
      possible_names <- c(paste0(ind_name, ":", shift_name), paste0(shift_name, ":", ind_name))
      int_name <- possible_names[possible_names %in% names(betas)][1]
      
      if (impact_model %in% c("full", "slope") && is.na(int_name)) stop(paste("Interaction term for disruption", k, "not found in model coefficients."))
      
      if (impact_model == "full") {
        var_lin <- cov_mat[ind_name, ind_name] + (b_k^2) * cov_mat[int_name, int_name] + (b_k * 2) * cov_mat[ind_name, int_name]
        lin_est <- betas[ind_name] + b_k * betas[int_name]
      } else if (impact_model == "level") {
        var_lin <- cov_mat[ind_name, ind_name]; lin_est <- betas[ind_name]
      } else if (impact_model == "slope") {
        var_lin <- (b_k^2) * cov_mat[int_name, int_name]; lin_est <- b_k * betas[int_name]
      }
      CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin)); TS <- lin_est / sqrt(var_lin)
      ret_vec[k, ] <- c(exp(lin_est), exp(CI_lin), round(2 * pnorm(abs(TS), lower.tail = FALSE), 2))
    }
    
  } else if (multi_model_type == "mixed") {
    fixed_eff <- lme4::fixef(model); rand_eff <- lme4::ranef(model)$disruption_id; cov_mat <- as.matrix(vcov(model))
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    
    if (impact_model == "full") {
      var_lin <- cov_mat["generic_indicator", "generic_indicator"] + (b_mult^2) * cov_mat["generic_shifted_time", "generic_shifted_time"] + (b_mult * 2) * cov_mat["generic_indicator", "generic_shifted_time"]
      lin_est <- fixed_eff["generic_indicator"] + b_mult * fixed_eff["generic_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat["generic_indicator", "generic_indicator"]; lin_est <- fixed_eff["generic_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat["generic_shifted_time", "generic_shifted_time"]; lin_est <- b_mult * fixed_eff["generic_shifted_time"]
    }
    
    RR_global <- exp(lin_est); CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin)); TS <- lin_est / sqrt(var_lin)
    ret_vec_global <- c(RR_global, exp(CI_lin), round(2 * pnorm(abs(TS), lower.tail = FALSE), 2))
    names(ret_vec_global) <- c("Global RR", "2.5% CI", "97.5% CI", "P-value")
    
    specific_mat <- matrix(NA, nrow = K, ncol = 1)
    rownames(specific_mat) <- paste0("Disruption_", 1:K)
    colnames(specific_mat) <- c("Specific RR (BLUP)")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      d_name <- paste0("Disruption_", k)
      if (d_name %in% rownames(rand_eff)) {
        if (impact_model == "full") blup <- (fixed_eff["generic_indicator"] + rand_eff[d_name, "generic_indicator"]) + b_k * (fixed_eff["generic_shifted_time"] + rand_eff[d_name, "generic_shifted_time"])
        else if (impact_model == "level") blup <- fixed_eff["generic_indicator"] + rand_eff[d_name, "generic_indicator"]
        else if (impact_model == "slope") blup <- b_k * (fixed_eff["generic_shifted_time"] + rand_eff[d_name, "generic_shifted_time"])
        specific_mat[k, 1] <- exp(blup)
      }
    }
    ret_vec <- list(Global_Average = ret_vec_global, Specific_Impacts = specific_mat)
  }
  
  # 4. OUTPUT & COUNTERFACTUALS 
  s <- summary(model)
  if (isTRUE(print_summary)) print(s)
  print(ret_vec)
  s[["RR"]] <- ret_vec 
  
  data$pred <- predict(model, type = "response")
  if (isTRUE(counterfactual)) {
    new_data <- data
    if (multi_model_type == "singular") {
      new_data$agg_indicator <- 0; new_data$agg_shifted_time <- 0
    } else if (multi_model_type == "multiple") {
      for (k in 1:K) { new_data[[paste0("indicator_", k)]] <- 0; new_data[[paste0("shifted_time_", k)]] <- 0 }
    } else if (multi_model_type == "mixed") {
      new_data$generic_indicator <- 0; new_data$generic_shifted_time <- 0
      new_data$disruption_id <- factor("Disruption_1", levels = levels(data$disruption_id))
    }
    
    data$predC <- predict(model, type = "response", newdata = new_data, allow.new.levels = TRUE)
    data$predC[!data$is_active] <- NA 
  }
  
  data <- tibble::as_tibble(data)
  return(list(model = model, model_summary = s, data = data))
}

#' ITS analysis for count outcomes with multiple interruptions and seasonal adjustment
#'
#' \code{its_poisson_fourier} fits a Poisson regression model adjusted for seasonality using Fourier terms. 
#' It accommodates multiple disruptions and returns the model, the summary (including the relative risk), 
#' and the original data appended with model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offset term.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL. Default value is NULL.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular" (homogenous effects) or "multiple" (distinct effects). Note: The "mixed" architecture is disabled here to prevent convergence failure with high-dimensional Fourier terms.
#' @param over_dispersion Logical - stating whether an over-dispersed quasi-Poisson model should be used (TRUE), or a regular Poisson (FALSE). Default is FALSE.
#' @param freq A positive integer describing the frequency of the time series (e.g., 12 for monthly data).
#' @param keep_significant_fourier Logical - indicating whether only the significant Fourier terms should be retained. Default is TRUE (the model is fitted twice to filter insignificant frequencies).
#' @param impact_model A string specifying the assumed impact model: "full", "level", or "slope". Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should be returned. Default is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed. Default is FALSE.
#' @return The function returns a list with three elements: the fitted Poisson regression model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' fit <- its_poisson_fourier(data=data, form=form, offset_name="labour", time_name="time", 
#'                            intervention_start_indices=start_indices, multi_model_type="multiple", 
#'                            over_dispersion=TRUE, freq=12, keep_significant_fourier=TRUE, 
#'                            impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export

its_poisson_fourier <- function(data, form, offset_name=NULL, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", over_dispersion=FALSE, freq, keep_significant_fourier=TRUE, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  pred <- predC <- indicator <- NULL
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(rlang::is_formula(form))) stop("Please make sure that form is a formula object")
  if (!(time_name %in% colnames(data))) stop("Please make sure that time_name belongs to colnames(data)")
  
  if (!all(data[[time_name]] == (1:nrow(data)))) {
    data[[time_name]] = 1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  
  if (!is.null(offset_name)) {
    if (!(offset_name %in% colnames(data))) stop("Please make sure that offset_name belongs to colnames(data)")
    if (offset_name %in% covariates | paste0("offset(log(",offset_name,"))") %in% covariates | paste0("offset(",offset_name,")") %in% covariates | paste0("log(",offset_name,")") %in% covariates) {
      stop("The offset term should not be included in the formula. Please supply a formula without the offset, and add the name of the offset column using the argument offset_name.")
    }
  }
  
  if (!is.logical(over_dispersion)) stop("over_dispersion must be either TRUE or FALSE")
  if (!(impact_model %in% c("full","level","slope"))) stop("Please enter a valid impact_model name. Possible impact models are \"full\", \"level\", and \"slope\".")
  if (!(response %in% colnames(data))) stop("Please make sure that the response variable on the left hand side of the formula object belongs to colnames(data)")
  if (!all(covariates %in% colnames(data))) stop("Please make sure that the covariates on the right hand side of the formula object belong to colnames(data)")
  if (!is.logical(keep_significant_fourier)) stop("keep_significant_fourier must be either TRUE or FALSE")
  
  n <- nrow(data)
  if (!(freq > 1 & freq <= n)) stop("Please make sure that the freq is a value greater than 1 and less or equal to nrow(data)")
  K <- length(intervention_start_indices)
  
  # 1. Helper & Dynamic Formula
  data <- build_disruption_features(data, time_name, intervention_start_indices, intervention_lengths)
  
  if (multi_model_type == "singular") {
    data$agg_indicator <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE])
    data$agg_shifted_time <- 0
    for (k in 1:K) data$agg_shifted_time <- data$agg_shifted_time + (data[[paste0("indicator_", k)]] * data[[paste0("shifted_time_", k)]])
    
    if (impact_model == "full") form_update <- update(form, as.formula("~ . + agg_indicator + agg_shifted_time"))
    else if (impact_model == "level") form_update <- update(form, as.formula("~ . + agg_indicator"))
    else if (impact_model == "slope") form_update <- update(form, as.formula("~ . + agg_shifted_time"))
    
  } else if (multi_model_type == "multiple") {
    terms <- c()
    for (k in 1:K) {
      ind <- paste0("indicator_", k)
      shift <- paste0("shifted_time_", k)
      if (impact_model == "full") terms <- c(terms, ind, paste0(ind, ":", shift))
      else if (impact_model == "level") terms <- c(terms, ind)
      else if (impact_model == "slope") terms <- c(terms, paste0(ind, ":", shift))
    }
    form_update <- update(form, as.formula(paste("~ . +", paste(terms, collapse = " + "))))
  } else stop("Please enter a valid multi_model_type: 'singular', 'multiple', or 'mixed'.")
  
  if (!is.null(offset_name)) form_update <- update(form_update, as.formula(paste0("~ . + offset(log(", sym(offset_name), "))")))
  
  # 2. Seasonality (Fourier) Logic
  ts_outcome <- ts(data[[response]], frequency = freq)
  fourier_mat <- forecast::fourier(ts_outcome, K = freq/2)
  colnames(fourier_mat) <- gsub("-", ".", colnames(fourier_mat))
  fourier_df <- data.frame(fourier_mat)
  colnames(fourier_df) <- colnames(fourier_mat)
  data <- cbind(data, fourier_df)
  form_update_full <- update(form_update, as.formula(paste0("~ . +", paste(colnames(fourier_mat), collapse= "+"))))
  
  set.seed(1)
  model <- glm(form_update_full, data, family = ifelse(!isTRUE(over_dispersion), poisson, quasipoisson))
  
  if (isTRUE(keep_significant_fourier)) {
    all_significant_terms <- summary(model)$coeff[,4] < 0.05
    
    # Dynamically build exclusion list so we don't drop intervention terms
    exclude_terms <- c("(Intercept)", time_name)
    if (multi_model_type == "singular") {
      exclude_terms <- c(exclude_terms, "agg_indicator", "agg_shifted_time")
    } else {
      for (k in 1:K) exclude_terms <- c(exclude_terms, paste0("indicator_", k), paste0("shifted_time_", k), paste0("indicator_", k, ":shifted_time_", k))
    }
    
    all_fourier_terms <- all_significant_terms[!names(all_significant_terms) %in% exclude_terms]
    selected_fourier <- names(all_fourier_terms)[all_fourier_terms == TRUE]
    form_update_sig <- update(form_update, as.formula(paste0("~ . +", paste(selected_fourier, collapse= "+"))))
    
    set.seed(1)
    model <- glm(form_update_sig, data, family = ifelse(!isTRUE(over_dispersion), poisson, quasipoisson))
  }
  
  # 3. Effect Size & Inference
  if (multi_model_type == "singular") data$is_active <- data$agg_indicator > 0
  else data$is_active <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE]) > 0
  
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  total_L <- sum(intervention_lengths)
  cov_mat <- vcov(model)
  
  if (multi_model_type == "singular") {
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    if (impact_model == "full") {
      var_lin <- cov_mat["agg_indicator", "agg_indicator"] + (b_mult^2) * cov_mat["agg_shifted_time", "agg_shifted_time"] + (b_mult * 2) * cov_mat["agg_indicator", "agg_shifted_time"]
      lin_est <- model$coefficients["agg_indicator"] + b_mult * model$coefficients["agg_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat["agg_indicator", "agg_indicator"]
      lin_est <- model$coefficients["agg_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat["agg_shifted_time", "agg_shifted_time"]
      lin_est <- b_mult * model$coefficients["agg_shifted_time"]
    }
    RR <- exp(lin_est)
    CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
    CI_RR <- exp(CI_lin)
    TS <- lin_est / sqrt(var_lin)
    p_value <- round(2 * pnorm(abs(TS), lower.tail = FALSE), 2)
    ret_vec <- c(RR, CI_RR, p_value)
    names(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
  } else if (multi_model_type == "multiple") {
    ret_vec <- matrix(NA, nrow = K, ncol = 4)
    rownames(ret_vec) <- paste0("Disruption_", 1:K)
    colnames(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      ind_name <- paste0("indicator_", k)
      shift_name <- paste0("shifted_time_", k)
      
      possible_names <- c(paste0(ind_name, ":", shift_name), paste0(shift_name, ":", ind_name))
      int_name <- possible_names[possible_names %in% names(model$coefficients)][1]
      
      if (impact_model %in% c("full", "slope") && is.na(int_name)) {
        stop(paste("Interaction term for disruption", k, "not found in model coefficients."))
      }
      
      if (impact_model == "full") {
        var_lin <- cov_mat[ind_name, ind_name] + (b_k^2) * cov_mat[int_name, int_name] + (b_k * 2) * cov_mat[ind_name, int_name]
        lin_est <- model$coefficients[ind_name] + b_k * model$coefficients[int_name]
      } else if (impact_model == "level") {
        var_lin <- cov_mat[ind_name, ind_name]
        lin_est <- model$coefficients[ind_name]
      } else if (impact_model == "slope") {
        var_lin <- (b_k^2) * cov_mat[int_name, int_name]
        lin_est <- b_k * model$coefficients[int_name]
      }
      RR_k <- exp(lin_est)
      CI_lin_k <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
      TS_k <- lin_est / sqrt(var_lin)
      ret_vec[k, ] <- c(RR_k, exp(CI_lin_k), round(2 * pnorm(abs(TS_k), lower.tail = FALSE), 2))
    }
  }
  
  s <- summary(model)
  if (isTRUE(print_summary)) print(s)
  print(ret_vec)
  s[["RR"]] <- ret_vec
  
  # 4. Predictions & Counterfactuals
  data$pred <- predict(model, type = "response")
  if (isTRUE(counterfactual)) {
    new_data <- data
    if (multi_model_type == "singular") {
      new_data$agg_indicator <- 0
      new_data$agg_shifted_time <- 0
    } else {
      for (k in 1:K) {
        new_data[[paste0("indicator_", k)]] <- 0
        new_data[[paste0("shifted_time_", k)]] <- 0
      }
    }
    data$predC <- predict(model, type = "response", newdata = new_data)
    data$predC[!data$is_active] <- NA 
  }
  return(list(model = model, model_summary = s, data = tibble::as_tibble(data)))
}

#' ITS analysis for zero-inflated count outcomes with multiple interruptions
#'
#' \code{its_zero_inflated} fits a Zero-Inflated Poisson (ZIP) regression model to an ITS with one or more disruption events. 
#' It models structural zeros via an intercept-only logistic component, and intervention effects via the count component.
#' Returns the model, the summary (including Relative Risk), and the original data appended with model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right. Should not contain an offset term.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL. Default value is NULL.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular", "multiple", or "mixed" (which relies on glmmTMB). Default is "singular".
#' @param impact_model A string specifying the assumed impact model: "full", "level", or "slope". Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned. Default is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed. Default is FALSE.
#' @return The function returns a list with three elements: the fitted zero-inflated Poisson model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- zero_inflated_sim_data
#' form <- as.formula("monthly_total ~ time")
#' start_indices <- c(which(data$Year==2020 & data$Month==3), which(data$Year==2021 & data$Month==1))
#' fit <- its_zero_inflated(data=data, form=form, time_name="time", 
#'                          intervention_start_indices=start_indices, 
#'                          multi_model_type="multiple", impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom pscl zeroinfl
#' @importFrom MASS mvrnorm
#' @export

its_zero_inflated <- function(data, form, offset_name=NULL, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  pred <- predC <- indicator <- NULL
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(rlang::is_formula(form))) stop("Please make sure that form is a formula object")
  if (!(time_name %in% colnames(data))) stop("Please make sure that time_name belongs to colnames(data)")
  
  if (!all(data[[time_name]] == (1:nrow(data)))) {
    data[[time_name]] = 1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(impact_model %in% c("full","level","slope"))) stop("Please enter a valid impact_model name: 'full', 'level', or 'slope'.")
  
  n <- nrow(data)
  K <- length(intervention_start_indices)
  
  # 1. DATA PREP & DYNAMIC FORMULA BUILDER 
  if (multi_model_type %in% c("singular", "multiple")) {
    data <- build_disruption_features(data, time_name, intervention_start_indices, intervention_lengths)
  } else if (multi_model_type == "mixed") {
    data <- build_mixed_features(data, time_name, intervention_start_indices, intervention_lengths)
  }
  
  # Base formula string extraction
  base_form_str <- paste(as.character(form)[2], "~", as.character(form)[3])
  if (!is.null(offset_name)) base_form_str <- paste0(base_form_str, " + offset(log(", offset_name, "))")
  
  if (multi_model_type == "singular") {
    data$agg_indicator <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE])
    data$agg_shifted_time <- 0
    for (k in 1:K) data$agg_shifted_time <- data$agg_shifted_time + (data[[paste0("indicator_", k)]] * data[[paste0("shifted_time_", k)]])
    
    if (impact_model == "full") form_str <- paste(base_form_str, "+ agg_indicator + agg_shifted_time | 1")
    else if (impact_model == "level") form_str <- paste(base_form_str, "+ agg_indicator | 1")
    else if (impact_model == "slope") form_str <- paste(base_form_str, "+ agg_shifted_time | 1")
    
  } else if (multi_model_type == "multiple") {
    terms <- c()
    for (k in 1:K) {
      ind <- paste0("indicator_", k)
      shift <- paste0("shifted_time_", k)
      if (impact_model == "full") terms <- c(terms, ind, paste0(ind, ":", shift))
      else if (impact_model == "level") terms <- c(terms, ind)
      else if (impact_model == "slope") terms <- c(terms, paste0(ind, ":", shift))
    }
    form_str <- paste(base_form_str, "+", paste(terms, collapse = " + "), "| 1")
    
  } else if (multi_model_type == "mixed") {
    if (impact_model == "full") form_str <- paste(base_form_str, "+ generic_indicator + generic_shifted_time + (0 + generic_indicator + generic_shifted_time | disruption_id)")
    else if (impact_model == "level") form_str <- paste(base_form_str, "+ generic_indicator + (0 + generic_indicator | disruption_id)")
    else if (impact_model == "slope") form_str <- paste(base_form_str, "+ generic_shifted_time + (0 + generic_shifted_time | disruption_id)")
  } else stop("Invalid multi_model_type: 'singular', 'multiple', or 'mixed'.")
  
  # 2. DUAL-ENGINE FITTING 
  set.seed(1)
  if (multi_model_type == "mixed") {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("The 'glmmTMB' package is required for the mixed zero-inflated model. Please install it.")
    model <- glmmTMB::glmmTMB(as.formula(form_str), data = data, ziformula = ~1, family = poisson)
  } else {
    model <- pscl::zeroinfl(as.formula(form_str), data = data, dist = "poisson")
  }
  
  # 3. DYNAMIC PARAMETER EXTRACTION 
  if (multi_model_type %in% c("singular", "multiple")) {
    if (multi_model_type == "singular") data$is_active <- data$agg_indicator > 0
    else data$is_active <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE]) > 0
  } else {
    data$is_active <- data$generic_indicator > 0
  }
  
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  total_L <- sum(intervention_lengths)
  
  # Route the coefficient & covariance extraction based on engine
  if (multi_model_type == "mixed") {
    betas <- glmmTMB::fixef(model)$cond
    cov_mat <- vcov(model)$cond
    prefix <- ""
  } else {
    betas <- model$coefficients$count
    cov_mat <- vcov(model)
    prefix <- "count_"
  }
  
  # 4. EFFECT SIZE (Relative Risk) 
  if (multi_model_type == "singular") {
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    if (impact_model == "full") {
      var_lin <- cov_mat[paste0(prefix, "agg_indicator"), paste0(prefix, "agg_indicator")] + 
        (b_mult^2) * cov_mat[paste0(prefix, "agg_shifted_time"), paste0(prefix, "agg_shifted_time")] + 
        (b_mult * 2) * cov_mat[paste0(prefix, "agg_indicator"), paste0(prefix, "agg_shifted_time")]
      lin_est <- betas["agg_indicator"] + b_mult * betas["agg_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat[paste0(prefix, "agg_indicator"), paste0(prefix, "agg_indicator")]
      lin_est <- betas["agg_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat[paste0(prefix, "agg_shifted_time"), paste0(prefix, "agg_shifted_time")]
      lin_est <- b_mult * betas["agg_shifted_time"]
    }
    
    RR <- exp(lin_est)
    CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
    TS <- lin_est / sqrt(var_lin)
    p_value <- round(2 * pnorm(abs(TS), lower.tail = FALSE), 2)
    ret_vec <- c(RR, exp(CI_lin), p_value)
    names(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
  } else if (multi_model_type == "multiple") {
    ret_vec <- matrix(NA, nrow = K, ncol = 4)
    rownames(ret_vec) <- paste0("Disruption_", 1:K)
    colnames(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      ind_name <- paste0("indicator_", k)
      shift_name <- paste0("shifted_time_", k)
      
      possible_names <- c(paste0(ind_name, ":", shift_name), paste0(shift_name, ":", ind_name))
      int_name <- possible_names[possible_names %in% names(betas)][1]
      
      if (impact_model %in% c("full", "slope") && is.na(int_name)) stop(paste("Interaction term for disruption", k, "not found in model coefficients."))
      
      if (impact_model == "full") {
        var_lin <- cov_mat[paste0(prefix, ind_name), paste0(prefix, ind_name)] + 
          (b_k^2) * cov_mat[paste0(prefix, int_name), paste0(prefix, int_name)] + 
          (b_k * 2) * cov_mat[paste0(prefix, ind_name), paste0(prefix, int_name)]
        lin_est <- betas[ind_name] + b_k * betas[int_name]
      } else if (impact_model == "level") {
        var_lin <- cov_mat[paste0(prefix, ind_name), paste0(prefix, ind_name)]
        lin_est <- betas[ind_name]
      } else if (impact_model == "slope") {
        var_lin <- (b_k^2) * cov_mat[paste0(prefix, int_name), paste0(prefix, int_name)]
        lin_est <- b_k * betas[int_name]
      }
      CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
      TS <- lin_est / sqrt(var_lin)
      ret_vec[k, ] <- c(exp(lin_est), exp(CI_lin), round(2 * pnorm(abs(TS), lower.tail = FALSE), 2))
    }
    
  } else if (multi_model_type == "mixed") {
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    if (impact_model == "full") {
      var_lin <- cov_mat["generic_indicator", "generic_indicator"] + (b_mult^2) * cov_mat["generic_shifted_time", "generic_shifted_time"] + (b_mult * 2) * cov_mat["generic_indicator", "generic_shifted_time"]
      lin_est <- betas["generic_indicator"] + b_mult * betas["generic_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat["generic_indicator", "generic_indicator"]
      lin_est <- betas["generic_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat["generic_shifted_time", "generic_shifted_time"]
      lin_est <- b_mult * betas["generic_shifted_time"]
    }
    
    RR_global <- exp(lin_est)
    CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
    TS <- lin_est / sqrt(var_lin)
    ret_vec_global <- c(RR_global, exp(CI_lin), round(2 * pnorm(abs(TS), lower.tail = FALSE), 2))
    names(ret_vec_global) <- c("Global RR", "2.5% CI", "97.5% CI", "P-value")
    
    rand_eff <- glmmTMB::ranef(model)$cond$disruption_id
    specific_mat <- matrix(NA, nrow = K, ncol = 1)
    rownames(specific_mat) <- paste0("Disruption_", 1:K)
    colnames(specific_mat) <- c("Specific RR (BLUP)")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      d_name <- paste0("Disruption_", k)
      if (d_name %in% rownames(rand_eff)) {
        if (impact_model == "full") blup <- (betas["generic_indicator"] + rand_eff[d_name, "generic_indicator"]) + b_k * (betas["generic_shifted_time"] + rand_eff[d_name, "generic_shifted_time"])
        else if (impact_model == "level") blup <- betas["generic_indicator"] + rand_eff[d_name, "generic_indicator"]
        else if (impact_model == "slope") blup <- b_k * (betas["generic_shifted_time"] + rand_eff[d_name, "generic_shifted_time"])
        specific_mat[k, 1] <- exp(blup)
      }
    }
    ret_vec <- list(Global_Average = ret_vec_global, Specific_Impacts = specific_mat)
  }
  
  # 5. OUTPUT & COUNTERFACTUALS 
  s <- summary(model)
  if (isTRUE(print_summary)) print(s)
  print(ret_vec)
  s[["RR"]] <- ret_vec 
  
  data$pred <- predict(model, type = "response")
  if (isTRUE(counterfactual)) {
    new_data <- data
    if (multi_model_type == "singular") {
      new_data$agg_indicator <- 0; new_data$agg_shifted_time <- 0
    } else if (multi_model_type == "multiple") {
      for (k in 1:K) { new_data[[paste0("indicator_", k)]] <- 0; new_data[[paste0("shifted_time_", k)]] <- 0 }
    } else if (multi_model_type == "mixed") {
      new_data$generic_indicator <- 0; new_data$generic_shifted_time <- 0
      # Matching the factor level generated by helper.R
      new_data$disruption_id <- factor("Disruption_1", levels = levels(data$disruption_id))
    }
    
    if (multi_model_type == "mixed") data$predC <- predict(model, type = "response", newdata = new_data, allow.new.levels = TRUE)
    else data$predC <- predict(model, type = "response", newdata = new_data)
    data$predC[!data$is_active] <- NA 
  }
  
  data <- tibble::as_tibble(data)
  return(list(model = model, model_summary = s, data = data))
}

#' ITS analysis for zero-inflated count outcomes with seasonal adjustment
#'
#' \code{its_zero_inflated_fourier} fits a zero-inflated Poisson regression model adjusted for seasonality 
#' to an ITS with multiple interruptions. It returns the fitted model, the summary (including relative risk), 
#' and the original data appended with model predictions.
#' Note: Combining mixed-effects multi-disruption architectures with Fourier terms and zero-inflation is 
#' computationally heavy and may be prone to non-convergence.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right. Should not contain an offset term.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL. Default value is NULL.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular", "multiple", or "mixed". Default is "singular".
#' @param freq A positive integer describing the frequency of the time series (e.g., 12 for monthly data).
#' @param impact_model A string specifying the assumed impact model: "full", "level", or "slope". Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned. Default is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed. Default is FALSE.
#' @return The function returns a list with three elements: the fitted zero-inflated Poisson model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- zero_inflated_sim_data
#' form <- as.formula("monthly_total ~ time")
#' start_indices <- c(which(data$Year==2020 & data$Month==3), which(data$Year==2021 & data$Month==1))
#' fit <- its_zero_inflated_fourier(data=data, form=form, time_name="time", 
#'                                  intervention_start_indices=start_indices, freq=12, 
#'                                  multi_model_type="multiple", impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom pscl zeroinfl
#' @importFrom MASS mvrnorm
#' @importFrom forecast fourier
#' @export

its_zero_inflated_fourier <- function(data, form, offset_name=NULL, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", freq, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  pred <- predC <- indicator <- NULL
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(rlang::is_formula(form))) stop("Please make sure that form is a formula object")
  if (!(time_name %in% colnames(data))) stop("Please make sure that time_name belongs to colnames(data)")
  
  if (!all(data[[time_name]] == (1:nrow(data)))) {
    data[[time_name]] = 1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(impact_model %in% c("full","level","slope"))) stop("Please enter a valid impact_model name: 'full', 'level', or 'slope'.")
  
  n <- nrow(data)
  K <- length(intervention_start_indices)
  if (!(freq > 1 & freq <= n)) stop("Please make sure that the freq is a value greater than 1 and less or equal to nrow(data)")
  
  # 1. DATA PREP & DYNAMIC FORMULA BUILDER 
  if (multi_model_type %in% c("singular", "multiple")) {
    data <- build_disruption_features(data, time_name, intervention_start_indices, intervention_lengths)
  } else if (multi_model_type == "mixed") {
    data <- build_mixed_features(data, time_name, intervention_start_indices, intervention_lengths)
  }
  
  # Generate and append Fourier terms
  ts_outcome <- ts(data[[response]], frequency = freq)
  fourier_mat <- forecast::fourier(ts_outcome, K = freq/2)
  colnames(fourier_mat) <- gsub("-", ".", colnames(fourier_mat))
  data <- cbind(data, data.frame(fourier_mat))
  fourier_str <- paste(colnames(fourier_mat), collapse = " + ")
  
  base_form_str <- paste(as.character(form)[2], "~", as.character(form)[3], "+", fourier_str)
  if (!is.null(offset_name)) base_form_str <- paste0(base_form_str, " + offset(log(", offset_name, "))")
  
  if (multi_model_type == "singular") {
    data$agg_indicator <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE])
    data$agg_shifted_time <- 0
    for (k in 1:K) data$agg_shifted_time <- data$agg_shifted_time + (data[[paste0("indicator_", k)]] * data[[paste0("shifted_time_", k)]])
    
    if (impact_model == "full") form_str <- paste(base_form_str, "+ agg_indicator + agg_shifted_time | 1")
    else if (impact_model == "level") form_str <- paste(base_form_str, "+ agg_indicator | 1")
    else if (impact_model == "slope") form_str <- paste(base_form_str, "+ agg_shifted_time | 1")
    
  } else if (multi_model_type == "multiple") {
    terms <- c()
    for (k in 1:K) {
      ind <- paste0("indicator_", k); shift <- paste0("shifted_time_", k)
      if (impact_model == "full") terms <- c(terms, ind, paste0(ind, ":", shift))
      else if (impact_model == "level") terms <- c(terms, ind)
      else if (impact_model == "slope") terms <- c(terms, paste0(ind, ":", shift))
    }
    form_str <- paste(base_form_str, "+", paste(terms, collapse = " + "), "| 1")
    
  } else if (multi_model_type == "mixed") {
    if (impact_model == "full") form_str <- paste(base_form_str, "+ generic_indicator + generic_shifted_time + (0 + generic_indicator + generic_shifted_time | disruption_id)")
    else if (impact_model == "level") form_str <- paste(base_form_str, "+ generic_indicator + (0 + generic_indicator | disruption_id)")
    else if (impact_model == "slope") form_str <- paste(base_form_str, "+ generic_shifted_time + (0 + generic_shifted_time | disruption_id)")
  } 
  
  # 2. DUAL-ENGINE FITTING 
  set.seed(1)
  if (multi_model_type == "mixed") {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("The 'glmmTMB' package is required for the mixed zero-inflated model. Please install it.")
    model <- glmmTMB::glmmTMB(as.formula(form_str), data = data, ziformula = ~1, family = poisson)
  } else {
    model <- pscl::zeroinfl(as.formula(form_str), data = data, dist = "poisson")
  }
  
  # 3. DYNAMIC PARAMETER EXTRACTION 
  if (multi_model_type %in% c("singular", "multiple")) {
    if (multi_model_type == "singular") data$is_active <- data$agg_indicator > 0
    else data$is_active <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE]) > 0
  } else {
    data$is_active <- data$generic_indicator > 0
  }
  
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  total_L <- sum(intervention_lengths)
  
  if (multi_model_type == "mixed") {
    betas <- glmmTMB::fixef(model)$cond
    cov_mat <- vcov(model)$cond
    prefix <- ""
  } else {
    betas <- model$coefficients$count
    cov_mat <- vcov(model)
    prefix <- "count_"
  }
  
  # 4. EFFECT SIZE (Relative Risk) 
  if (multi_model_type == "singular") {
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    if (impact_model == "full") {
      var_lin <- cov_mat[paste0(prefix, "agg_indicator"), paste0(prefix, "agg_indicator")] + (b_mult^2) * cov_mat[paste0(prefix, "agg_shifted_time"), paste0(prefix, "agg_shifted_time")] + (b_mult * 2) * cov_mat[paste0(prefix, "agg_indicator"), paste0(prefix, "agg_shifted_time")]
      lin_est <- betas["agg_indicator"] + b_mult * betas["agg_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat[paste0(prefix, "agg_indicator"), paste0(prefix, "agg_indicator")]
      lin_est <- betas["agg_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat[paste0(prefix, "agg_shifted_time"), paste0(prefix, "agg_shifted_time")]
      lin_est <- b_mult * betas["agg_shifted_time"]
    }
    
    RR <- exp(lin_est); CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
    TS <- lin_est / sqrt(var_lin); p_value <- round(2 * pnorm(abs(TS), lower.tail = FALSE), 2)
    ret_vec <- c(RR, exp(CI_lin), p_value); names(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
  } else if (multi_model_type == "multiple") {
    ret_vec <- matrix(NA, nrow = K, ncol = 4)
    rownames(ret_vec) <- paste0("Disruption_", 1:K)
    colnames(ret_vec) <- c("RR", "2.5% CI", "97.5% CI", "P-value")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      ind_name <- paste0("indicator_", k)
      shift_name <- paste0("shifted_time_", k)
      
      possible_names <- c(paste0(ind_name, ":", shift_name), paste0(shift_name, ":", ind_name))
      int_name <- possible_names[possible_names %in% names(betas)][1]
      
      if (impact_model %in% c("full", "slope") && is.na(int_name)) stop(paste("Interaction term for disruption", k, "not found in model coefficients."))
      
      if (impact_model == "full") {
        var_lin <- cov_mat[paste0(prefix, ind_name), paste0(prefix, ind_name)] + 
          (b_k^2) * cov_mat[paste0(prefix, int_name), paste0(prefix, int_name)] + 
          (b_k * 2) * cov_mat[paste0(prefix, ind_name), paste0(prefix, int_name)]
        lin_est <- betas[ind_name] + b_k * betas[int_name]
      } else if (impact_model == "level") {
        var_lin <- cov_mat[paste0(prefix, ind_name), paste0(prefix, ind_name)]
        lin_est <- betas[ind_name]
      } else if (impact_model == "slope") {
        var_lin <- (b_k^2) * cov_mat[paste0(prefix, int_name), paste0(prefix, int_name)]
        lin_est <- b_k * betas[int_name]
      }
      CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin))
      TS <- lin_est / sqrt(var_lin)
      ret_vec[k, ] <- c(exp(lin_est), exp(CI_lin), round(2 * pnorm(abs(TS), lower.tail = FALSE), 2))
    }
    
  } else if (multi_model_type == "mixed") {
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    if (impact_model == "full") {
      var_lin <- cov_mat["generic_indicator", "generic_indicator"] + (b_mult^2) * cov_mat["generic_shifted_time", "generic_shifted_time"] + (b_mult * 2) * cov_mat["generic_indicator", "generic_shifted_time"]
      lin_est <- betas["generic_indicator"] + b_mult * betas["generic_shifted_time"]
    } else if (impact_model == "level") {
      var_lin <- cov_mat["generic_indicator", "generic_indicator"]; lin_est <- betas["generic_indicator"]
    } else if (impact_model == "slope") {
      var_lin <- (b_mult^2) * cov_mat["generic_shifted_time", "generic_shifted_time"]; lin_est <- b_mult * betas["generic_shifted_time"]
    }
    
    RR_global <- exp(lin_est); CI_lin <- c(lin_est - 1.96 * sqrt(var_lin), lin_est + 1.96 * sqrt(var_lin)); TS <- lin_est / sqrt(var_lin)
    ret_vec_global <- c(RR_global, exp(CI_lin), round(2 * pnorm(abs(TS), lower.tail = FALSE), 2))
    names(ret_vec_global) <- c("Global RR", "2.5% CI", "97.5% CI", "P-value")
    
    rand_eff <- glmmTMB::ranef(model)$cond$disruption_id
    specific_mat <- matrix(NA, nrow = K, ncol = 1)
    rownames(specific_mat) <- paste0("Disruption_", 1:K)
    colnames(specific_mat) <- c("Specific RR (BLUP)")
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      d_name <- paste0("Disruption_", k)
      if (d_name %in% rownames(rand_eff)) {
        if (impact_model == "full") blup <- (betas["generic_indicator"] + rand_eff[d_name, "generic_indicator"]) + b_k * (betas["generic_shifted_time"] + rand_eff[d_name, "generic_shifted_time"])
        else if (impact_model == "level") blup <- betas["generic_indicator"] + rand_eff[d_name, "generic_indicator"]
        else if (impact_model == "slope") blup <- b_k * (betas["generic_shifted_time"] + rand_eff[d_name, "generic_shifted_time"])
        specific_mat[k, 1] <- exp(blup)
      }
    }
    ret_vec <- list(Global_Average = ret_vec_global, Specific_Impacts = specific_mat)
  }
  
  # 5. OUTPUT & COUNTERFACTUALS 
  s <- summary(model)
  if (isTRUE(print_summary)) print(s)
  print(ret_vec)
  s[["RR"]] <- ret_vec 
  
  data$pred <- predict(model, type = "response")
  if (isTRUE(counterfactual)) {
    new_data <- data
    if (multi_model_type == "singular") {
      new_data$agg_indicator <- 0; new_data$agg_shifted_time <- 0
    } else if (multi_model_type == "multiple") {
      for (k in 1:K) { new_data[[paste0("indicator_", k)]] <- 0; new_data[[paste0("shifted_time_", k)]] <- 0 }
    } else if (multi_model_type == "mixed") {
      new_data$generic_indicator <- 0; new_data$generic_shifted_time <- 0
      new_data$disruption_id <- factor("Disruption_1", levels = levels(data$disruption_id))
    }
    
    if (multi_model_type == "mixed") data$predC <- predict(model, type = "response", newdata = new_data, allow.new.levels = TRUE)
    else data$predC <- predict(model, type = "response", newdata = new_data)
    data$predC[!data$is_active] <- NA 
  }
  
  data <- tibble::as_tibble(data)
  return(list(model = model, model_summary = s, data = data))
}

#' ITS analysis for continuous outcomes with multiple interruptions (no seasonality)
#' 
#' \code{its_lm_wo_seas} fits a linear regression model to an ITS with one or more disruption events. 
#' It estimates the effect size via Cohen's d using a pooled standard deviation of factual and counterfactual predictions.
#' Confidence intervals and P-values are derived via parametric bootstrapping.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the continuous outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right. Should not contain an offset term.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent. Default is NULL.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular", "multiple", or "mixed". Default is "singular".
#' @param impact_model A string specifying the assumed impact model: "full", "level", or "slope". Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned in the data. Default is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed. Default is FALSE.
#' @return The function returns a list with three elements: the fitted linear regression model, the summary of the model (including Cohen's d), and the original data appended with model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' lengths <- c(6, 12)
#' fit <- its_lm_wo_seas(data=data, form=form, time_name="time", 
#'                       intervention_start_indices=start_indices, intervention_lengths=lengths, 
#'                       multi_model_type="multiple", impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm ts vcov update predict var na.omit quantile
#' @importFrom lme4 lmer fixef ranef bootMer
#' @importFrom MASS mvrnorm
#' @export

its_lm_wo_seas <- function(data, form, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  pred <- predC <- indicator <- NULL
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!rlang::is_formula(form)) stop("Please make sure that form is a formula object")
  if (!(time_name %in% colnames(data))) stop("Please make sure that time_name belongs to colnames(data)")
  
  if (!all(data[[time_name]] == (1:nrow(data)))) {
    data[[time_name]] = 1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(impact_model %in% c("full","level","slope"))) stop("Please enter a valid impact_model name: 'full', 'level', or 'slope'.")
  
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  
  n <- nrow(data)
  K <- length(intervention_start_indices)
  
  # 1. DATA PREP & FORMULA BUILDER 
  if (multi_model_type %in% c("singular", "multiple")) data <- build_disruption_features(data, time_name, intervention_start_indices, intervention_lengths)
  else if (multi_model_type == "mixed") data <- build_mixed_features(data, time_name, intervention_start_indices, intervention_lengths)
  
  if (multi_model_type == "singular") {
    data$agg_indicator <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE])
    data$agg_shifted_time <- 0
    for (k in 1:K) data$agg_shifted_time <- data$agg_shifted_time + (data[[paste0("indicator_", k)]] * data[[paste0("shifted_time_", k)]])
    
    if (impact_model == "full") form_update <- update(form, as.formula("~ . + agg_indicator + agg_shifted_time"))
    else if (impact_model == "level") form_update <- update(form, as.formula("~ . + agg_indicator"))
    else if (impact_model == "slope") form_update <- update(form, as.formula("~ . + agg_shifted_time"))
    
  } else if (multi_model_type == "multiple") {
    terms <- c()
    for (k in 1:K) {
      ind <- paste0("indicator_", k); shift <- paste0("shifted_time_", k)
      if (impact_model == "full") terms <- c(terms, ind, paste0(ind, ":", shift))
      else if (impact_model == "level") terms <- c(terms, ind)
      else if (impact_model == "slope") terms <- c(terms, paste0(ind, ":", shift))
    }
    form_update <- update(form, as.formula(paste("~ . +", paste(terms, collapse = " + "))))
    
  } else if (multi_model_type == "mixed") {
    if (impact_model == "full") form_update <- update(form, as.formula("~ . + generic_indicator + generic_shifted_time + (0 + generic_indicator + generic_shifted_time | disruption_id)"))
    else if (impact_model == "level") form_update <- update(form, as.formula("~ . + generic_indicator + (0 + generic_indicator | disruption_id)"))
    else if (impact_model == "slope") form_update <- update(form, as.formula("~ . + generic_shifted_time + (0 + generic_shifted_time | disruption_id)"))
  } 
  
  # 2. MODEL FITTING 
  if (multi_model_type == "mixed") {
    if (!requireNamespace("lme4", quietly = TRUE)) stop("The 'lme4' package is required for the mixed model.")
    model <- lme4::lmer(form_update, data)
  } else {
    model <- lm(form_update, data)
  }
  
  # 3. BASE PREDICTIONS (For Sp Calculation) 
  if (multi_model_type %in% c("singular", "multiple")) {
    if (multi_model_type == "singular") data$is_active <- data$agg_indicator > 0
    else data$is_active <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE]) > 0
  } else {
    data$is_active <- data$generic_indicator > 0
  }
  
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  total_L <- sum(intervention_lengths)
  
  data$pred <- predict(model, allow.new.levels = TRUE)
  new_data <- data
  if (multi_model_type == "singular") {
    new_data$agg_indicator <- 0; new_data$agg_shifted_time <- 0
  } else if (multi_model_type == "multiple") {
    for (k in 1:K) { new_data[[paste0("indicator_", k)]] <- 0; new_data[[paste0("shifted_time_", k)]] <- 0 }
  } else if (multi_model_type == "mixed") {
    new_data$generic_indicator <- 0; new_data$generic_shifted_time <- 0
    new_data$disruption_id <- factor("Disruption_1", levels = levels(data$disruption_id))
  }
  
  data$predC <- predict(model, newdata = new_data, allow.new.levels = TRUE)
  data$predC[!data$is_active] <- NA
  s1 <- var(data$pred[data$is_active]); s2 <- var(na.omit(data$predC)); Sp <- sqrt((s1 + s2) / 2)
  
  # 4. EFFECT SIZE (Cohen's d) & BOOTSTRAP INFERENCE 
  n_boot_samples <- 2000
  
  if (multi_model_type == "singular") {
    betas <- model$coefficients; cov_mat <- vcov(model); model_copy <- model
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    
    if (impact_model == "full") mean_diff <- betas["agg_indicator"] + b_mult * betas["agg_shifted_time"]
    else if (impact_model == "level") mean_diff <- betas["agg_indicator"]
    else if (impact_model == "slope") mean_diff <- b_mult * betas["agg_shifted_time"]
    
    cohen_d <- mean_diff / Sp
    cohen_d_mat <- matrix(NA, n_boot_samples, 2)
    
    for (i in 1:n_boot_samples) {
      set.seed(i)
      sampled_beta <- MASS::mvrnorm(n = 1, mu = betas, Sigma = cov_mat)
      sampled_beta_null <- sampled_beta
      
      if (impact_model == "full") {
        sampled_beta_null[c("agg_indicator", "agg_shifted_time")] <- sampled_beta[c("agg_indicator", "agg_shifted_time")] - betas[c("agg_indicator", "agg_shifted_time")]
        MD <- sampled_beta["agg_indicator"] + b_mult * sampled_beta["agg_shifted_time"]
        MD_null <- sampled_beta_null["agg_indicator"] + b_mult * sampled_beta_null["agg_shifted_time"]
      } else if (impact_model == "level") {
        sampled_beta_null["agg_indicator"] <- sampled_beta["agg_indicator"] - betas["agg_indicator"]; MD <- sampled_beta["agg_indicator"]; MD_null <- sampled_beta_null["agg_indicator"]
      } else if (impact_model == "slope") {
        sampled_beta_null["agg_shifted_time"] <- sampled_beta["agg_shifted_time"] - betas["agg_shifted_time"]; MD <- b_mult * sampled_beta["agg_shifted_time"]; MD_null <- b_mult * sampled_beta_null["agg_shifted_time"]
      }
      
      model_copy$coefficients <- sampled_beta; Sp_boot <- sqrt((var(predict(model_copy)[data$is_active]) + var(predict(model_copy, newdata=new_data)[data$is_active])) / 2)
      cohen_d_mat[i, 1] <- MD / Sp_boot
      
      model_copy$coefficients <- sampled_beta_null; Sp_boot_null <- sqrt((var(predict(model_copy)[data$is_active]) + var(predict(model_copy, newdata=new_data)[data$is_active])) / 2)
      cohen_d_mat[i, 2] <- MD_null / Sp_boot_null
    }
    
    ret_vec <- c(cohen_d, quantile(cohen_d_mat[, 1], 0.025), quantile(cohen_d_mat[, 1], 0.975), mean(abs(cohen_d_mat[, 2]) > abs(cohen_d)))
    names(ret_vec) <- c("Cohen's d", "2.5% CI", "97.5% CI", "P-value")
    
  } else if (multi_model_type == "multiple") {
    betas <- model$coefficients; cov_mat <- vcov(model); model_copy <- model
    ret_vec <- matrix(NA, nrow = K, ncol = 4)
    rownames(ret_vec) <- paste0("Disruption_", 1:K); colnames(ret_vec) <- c("Cohen's d", "2.5% CI", "97.5% CI", "P-value")
    
    cohen_boot_array <- array(NA, dim = c(n_boot_samples, 2, K)); base_cohen_d <- numeric(K)
    int_names_k <- character(K) # Store names to avoid searching 2000 times
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      ind_name <- paste0("indicator_", k); shift_name <- paste0("shifted_time_", k)
      
      # FIX: Defensive Interaction Mapping
      possible_names <- c(paste0(ind_name, ":", shift_name), paste0(shift_name, ":", ind_name))
      int_name <- possible_names[possible_names %in% names(betas)][1]
      int_names_k[k] <- int_name
      
      if (impact_model %in% c("full", "slope") && is.na(int_name)) stop(paste("Interaction term for disruption", k, "not found in model coefficients."))
      
      if (impact_model == "full") MD_k <- betas[ind_name] + b_k * betas[int_name]
      else if (impact_model == "level") MD_k <- betas[ind_name]
      else if (impact_model == "slope") MD_k <- b_k * betas[int_name]
      base_cohen_d[k] <- MD_k / Sp
    }
    
    for (i in 1:n_boot_samples) {
      set.seed(i)
      sampled_beta <- MASS::mvrnorm(n = 1, mu = betas, Sigma = cov_mat); sampled_beta_null <- sampled_beta
      
      for (k in 1:K) {
        int_name <- int_names_k[k]
        if (impact_model %in% c("full", "level")) sampled_beta_null[paste0("indicator_", k)] <- sampled_beta[paste0("indicator_", k)] - betas[paste0("indicator_", k)]
        if (impact_model %in% c("full", "slope")) sampled_beta_null[int_name] <- sampled_beta[int_name] - betas[int_name]
      }
      
      model_copy$coefficients <- sampled_beta; Sp_boot <- sqrt((var(predict(model_copy)[data$is_active]) + var(predict(model_copy, newdata=new_data)[data$is_active])) / 2)
      model_copy$coefficients <- sampled_beta_null; Sp_boot_null <- sqrt((var(predict(model_copy)[data$is_active]) + var(predict(model_copy, newdata=new_data)[data$is_active])) / 2)
      
      for (k in 1:K) {
        b_k <- (intervention_lengths[k] - 1) / 2
        ind_name <- paste0("indicator_", k); int_name <- int_names_k[k]
        
        if (impact_model == "full") {
          MD_k <- sampled_beta[ind_name] + b_k * sampled_beta[int_name]
          MD_null_k <- sampled_beta_null[ind_name] + b_k * sampled_beta_null[int_name]
        } else if (impact_model == "level") {
          MD_k <- sampled_beta[ind_name]; MD_null_k <- sampled_beta_null[ind_name]
        } else if (impact_model == "slope") {
          MD_k <- b_k * sampled_beta[int_name]; MD_null_k <- b_k * sampled_beta_null[int_name]
        }
        cohen_boot_array[i, 1, k] <- MD_k / Sp_boot; cohen_boot_array[i, 2, k] <- MD_null_k / Sp_boot_null
      }
    }
    for (k in 1:K) ret_vec[k, ] <- c(base_cohen_d[k], quantile(cohen_boot_array[, 1, k], 0.025), quantile(cohen_boot_array[, 1, k], 0.975), mean(abs(cohen_boot_array[, 2, k]) > abs(base_cohen_d[k])))
    
  } else if (multi_model_type == "mixed") {
    fixed_eff <- lme4::fixef(model); rand_eff <- lme4::ranef(model)$disruption_id
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    
    if (impact_model == "full") mean_diff_global <- fixed_eff["generic_indicator"] + b_mult * fixed_eff["generic_shifted_time"]
    else if (impact_model == "level") mean_diff_global <- fixed_eff["generic_indicator"]
    else if (impact_model == "slope") mean_diff_global <- b_mult * fixed_eff["generic_shifted_time"]
    cohen_d_global <- mean_diff_global / Sp
    
    boot_func <- function(fit) {
      f_eff <- lme4::fixef(fit); r_eff <- lme4::ranef(fit)$disruption_id; res <- numeric(1 + K) 
      if (impact_model == "full") res[1] <- f_eff["generic_indicator"] + b_mult * f_eff["generic_shifted_time"]
      else if (impact_model == "level") res[1] <- f_eff["generic_indicator"]
      else if (impact_model == "slope") res[1] <- b_mult * f_eff["generic_shifted_time"]
      res[1] <- res[1] / Sp
      
      for (k in 1:K) {
        b_k <- (intervention_lengths[k] - 1) / 2; d_name <- paste0("Disruption_", k)
        if (d_name %in% rownames(r_eff)) {
          if (impact_model == "full") blup <- (f_eff["generic_indicator"] + r_eff[d_name, "generic_indicator"]) + b_k * (f_eff["generic_shifted_time"] + r_eff[d_name, "generic_shifted_time"])
          else if (impact_model == "level") blup <- f_eff["generic_indicator"] + r_eff[d_name, "generic_indicator"]
          else if (impact_model == "slope") blup <- b_k * (f_eff["generic_shifted_time"] + r_eff[d_name, "generic_shifted_time"])
          res[1 + k] <- blup / Sp
        } else res[1 + k] <- NA
      }
      return(res)
    }
    
    boot_res <- lme4::bootMer(model, boot_func, nsim = n_boot_samples, type = "parametric")
    CI_low_global <- quantile(boot_res$t[, 1], 0.025, na.rm = TRUE); CI_up_global <- quantile(boot_res$t[, 1], 0.975, na.rm = TRUE)
    se_global <- sd(boot_res$t[, 1], na.rm = TRUE); p_value_global <- round(2 * pnorm(abs(cohen_d_global / se_global), lower.tail = FALSE), 2)
    ret_vec_global <- c(cohen_d_global, CI_low_global, CI_up_global, p_value_global); names(ret_vec_global) <- c("Global Cohen's d", "2.5% CI", "97.5% CI", "P-value")
    
    specific_mat <- matrix(NA, nrow = K, ncol = 4); rownames(specific_mat) <- paste0("Disruption_", 1:K); colnames(specific_mat) <- c("Specific Cohen's d (BLUP)", "2.5% CI", "97.5% CI", "P-value")
    for (k in 1:K) {
      d_k <- boot_res$t0[1 + k]; CI_low_k <- quantile(boot_res$t[, 1 + k], 0.025, na.rm = TRUE); CI_up_k <- quantile(boot_res$t[, 1 + k], 0.975, na.rm = TRUE)
      se_k <- sd(boot_res$t[, 1 + k], na.rm = TRUE); p_val_k <- round(2 * pnorm(abs(d_k / se_k), lower.tail = FALSE), 2)
      specific_mat[k, ] <- c(d_k, CI_low_k, CI_up_k, p_val_k)
    }
    ret_vec <- list(Global_Average = ret_vec_global, Specific_Impacts = specific_mat)
  }
  
  # 5. OUTPUT PROCESSING 
  summar <- summary(model)
  if (isTRUE(print_summary)) print(summar)
  print(ret_vec)
  summar[["Cohen's d"]] <- ret_vec 
  
  if (!isTRUE(counterfactual)) data$predC <- NULL
  data <- tibble::as_tibble(data)
  return(list(model = model, model_summary = summar, data = data))
}

#' ITS analysis for continuous outcomes with seasonal adjustments via Fourier terms
#'
#' \code{its_lm_fourier} fits a linear regression model adjusted for seasonality to an ITS with multiple interruptions. 
#' It estimates the effect size via Cohen's d using parametric bootstrapping. 
#' Note: The "mixed" multi-model architecture is disabled here to prevent convergence failure with high-dimensional Fourier terms.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the continuous outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right. Should not contain an offset term.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent. Default is NULL.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular" or "multiple". Default is "singular".
#' @param freq A positive integer describing the frequency of the time series (e.g., 12 for monthly data).
#' @param keep_significant_fourier Logical - indicating whether only the significant Fourier terms should be considered. Default is TRUE.
#' @param impact_model A string specifying the assumed impact model: "full", "level", or "slope". Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned. Default is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed. Default is FALSE.
#' @return The function returns a list with three elements: the fitted linear regression model, the summary of the model (including Cohen's d), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' fit <- its_lm_fourier(data=data, form=form, time_name="time", 
#'                       intervention_start_indices=start_indices, freq=12, 
#'                       multi_model_type="multiple", keep_significant_fourier=TRUE, 
#'                       impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export

its_lm_fourier <- function(data, form, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", freq, keep_significant_fourier=TRUE, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  pred <- predC <- indicator <- NULL
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(rlang::is_formula(form))) stop("Please make sure that form is a formula object")
  if (!(time_name %in% colnames(data))) stop("Please make sure that time_name belongs to colnames(data)")
  
  if (!all(data[[time_name]] == (1:nrow(data)))) {
    data[[time_name]] = 1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(impact_model %in% c("full","level","slope"))) stop("Please enter a valid impact_model name. Possible impact models are \"full\", \"level\", and \"slope\".")
  
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  if (!(response %in% colnames(data))) stop("Please make sure that the response variable on the left hand side of the formula object belongs to colnames(data)")
  if (!all(covariates %in% colnames(data))) stop("Please make sure that the covariates on the right hand side of the formula object belong to colnames(data)")
  
  n <- nrow(data)
  if (!(freq > 1 & freq <= n)) stop("Please make sure that the freq is a value greater than 1 and less or equal to nrow(data)")
  K <- length(intervention_start_indices)
  
  # 1. Helper & Dynamic Formula
  data <- build_disruption_features(data, time_name, intervention_start_indices, intervention_lengths)
  
  if (multi_model_type == "singular") {
    data$agg_indicator <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE])
    data$agg_shifted_time <- 0
    for (k in 1:K) data$agg_shifted_time <- data$agg_shifted_time + (data[[paste0("indicator_", k)]] * data[[paste0("shifted_time_", k)]])
    
    if (impact_model == "full") form_update <- update(form, as.formula("~ . + agg_indicator + agg_shifted_time"))
    else if (impact_model == "level") form_update <- update(form, as.formula("~ . + agg_indicator"))
    else if (impact_model == "slope") form_update <- update(form, as.formula("~ . + agg_shifted_time"))
    
  } else if (multi_model_type == "multiple") {
    terms <- c()
    for (k in 1:K) {
      ind <- paste0("indicator_", k)
      shift <- paste0("shifted_time_", k)
      if (impact_model == "full") terms <- c(terms, ind, paste0(ind, ":", shift))
      else if (impact_model == "level") terms <- c(terms, ind)
      else if (impact_model == "slope") terms <- c(terms, paste0(ind, ":", shift))
    }
    form_update <- update(form, as.formula(paste("~ . +", paste(terms, collapse = " + "))))
  } else stop("Please enter a valid multi_model_type: 'singular', 'multiple', or 'mixed'.")
  
  # 2. Seasonality (Fourier) Logic
  ts_outcome <- ts(data[[response]], frequency = freq)
  fourier_mat <- forecast::fourier(ts_outcome, K = freq/2)
  colnames(fourier_mat) <- gsub("-", ".", colnames(fourier_mat))
  fourier_df <- data.frame(fourier_mat)
  colnames(fourier_df) <- colnames(fourier_mat)
  data <- cbind(data, fourier_df)
  form_update_full <- update(form_update, as.formula(paste0("~ . +", paste(colnames(fourier_mat), collapse= "+"))))
  
  model <- lm(form_update_full, data)
  
  if (isTRUE(keep_significant_fourier)) {
    all_significant_terms <- summary(model)$coeff[,4] < 0.05
    
    exclude_terms <- c("(Intercept)", time_name)
    if (multi_model_type == "singular") {
      exclude_terms <- c(exclude_terms, "agg_indicator", "agg_shifted_time")
    } else {
      for (k in 1:K) exclude_terms <- c(exclude_terms, paste0("indicator_", k), paste0("shifted_time_", k), paste0("indicator_", k, ":shifted_time_", k))
    }
    
    all_fourier_terms <- all_significant_terms[!names(all_significant_terms) %in% exclude_terms]
    selected_fourier <- names(all_fourier_terms)[all_fourier_terms == TRUE]
    form_update_sig <- update(form_update, as.formula(paste0("~ . +", paste(selected_fourier, collapse= "+"))))
    model <- lm(form_update_sig, data)
  }
  
  # 3. Effect Size & Bootstrapping
  cov_mat <- vcov(model)
  if (multi_model_type == "singular") data$is_active <- data$agg_indicator > 0
  else data$is_active <- rowSums(data[, paste0("indicator_", 1:K), drop = FALSE]) > 0
  
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  total_L <- sum(intervention_lengths)
  
  data$pred <- predict(model, type="response")
  new_data <- data
  if (multi_model_type == "singular") {
    new_data$agg_indicator <- 0
    new_data$agg_shifted_time <- 0
  } else {
    for (k in 1:K) {
      new_data[[paste0("indicator_", k)]] <- 0
      new_data[[paste0("shifted_time_", k)]] <- 0
    }
  }
  data$predC <- predict(model, type="response", newdata=new_data)
  data$predC[!data$is_active] <- NA
  
  s1 <- var(data$pred[data$is_active])
  s2 <- var(na.omit(data$predC))
  Sp <- sqrt((s1 + s2) / 2)
  
  betas <- model$coefficients
  n_boot_samples <- 2000
  model_copy <- model
  
  if (multi_model_type == "singular") {
    b_mult <- sum((intervention_lengths * (intervention_lengths - 1)) / 2) / total_L
    if (impact_model == "full") mean_diff <- betas["agg_indicator"] + b_mult * betas["agg_shifted_time"]
    else if (impact_model == "level") mean_diff <- betas["agg_indicator"]
    else if (impact_model == "slope") mean_diff <- b_mult * betas["agg_shifted_time"]
    
    cohen_d <- mean_diff / Sp
    cohen_d_mat <- matrix(NA, n_boot_samples, 2)
    
    for (i in 1:n_boot_samples) {
      set.seed(i)
      sampled_beta <- MASS::mvrnorm(n = 1, mu = betas, Sigma = cov_mat)
      sampled_beta_null <- sampled_beta
      
      if (impact_model == "full") {
        sampled_beta_null[c("agg_indicator", "agg_shifted_time")] <- sampled_beta[c("agg_indicator", "agg_shifted_time")] - betas[c("agg_indicator", "agg_shifted_time")]
        MD <- sampled_beta["agg_indicator"] + b_mult * sampled_beta["agg_shifted_time"]
        MD_null <- sampled_beta_null["agg_indicator"] + b_mult * sampled_beta_null["agg_shifted_time"]
      } else if (impact_model == "level") {
        sampled_beta_null["agg_indicator"] <- sampled_beta["agg_indicator"] - betas["agg_indicator"]
        MD <- sampled_beta["agg_indicator"]
        MD_null <- sampled_beta_null["agg_indicator"]
      } else if (impact_model == "slope") {
        sampled_beta_null["agg_shifted_time"] <- sampled_beta["agg_shifted_time"] - betas["agg_shifted_time"]
        MD <- b_mult * sampled_beta["agg_shifted_time"]
        MD_null <- b_mult * sampled_beta_null["agg_shifted_time"]
      }
      
      model_copy$coefficients <- sampled_beta
      Sp_boot <- sqrt((var(predict(model_copy, type="response")[data$is_active]) + var(predict(model_copy, type="response", newdata=new_data)[data$is_active])) / 2)
      cohen_d_mat[i, 1] <- MD / Sp_boot
      
      model_copy$coefficients <- sampled_beta_null
      Sp_boot_null <- sqrt((var(predict(model_copy, type="response")[data$is_active]) + var(predict(model_copy, type="response", newdata=new_data)[data$is_active])) / 2)
      cohen_d_mat[i, 2] <- MD_null / Sp_boot_null
    }
    
    ret_vec <- c(cohen_d, quantile(cohen_d_mat[, 1], 0.025), quantile(cohen_d_mat[, 1], 0.975), mean(abs(cohen_d_mat[, 2]) > abs(cohen_d)))
    names(ret_vec) <- c("Cohen's d", "2.5% CI", "97.5% CI", "P-value")
    
  } else if (multi_model_type == "multiple") {
    ret_vec <- matrix(NA, nrow = K, ncol = 4)
    rownames(ret_vec) <- paste0("Disruption_", 1:K)
    colnames(ret_vec) <- c("Cohen's d", "2.5% CI", "97.5% CI", "P-value")
    cohen_boot_array <- array(NA, dim = c(n_boot_samples, 2, K))
    base_cohen_d <- numeric(K)
    int_names_k <- character(K) # Store names to avoid searching 2000 times
    
    for (k in 1:K) {
      b_k <- (intervention_lengths[k] - 1) / 2
      ind_name <- paste0("indicator_", k)
      shift_name <- paste0("shifted_time_", k)
      
      # VITAL FIX: Defensive Interaction Mapping
      possible_names <- c(paste0(ind_name, ":", shift_name), paste0(shift_name, ":", ind_name))
      int_name <- possible_names[possible_names %in% names(betas)][1]
      int_names_k[k] <- int_name
      
      if (impact_model %in% c("full", "slope") && is.na(int_name)) stop(paste("Interaction term for disruption", k, "not found in model coefficients."))
      
      if (impact_model == "full") MD_k <- betas[ind_name] + b_k * betas[int_name]
      else if (impact_model == "level") MD_k <- betas[ind_name]
      else if (impact_model == "slope") MD_k <- b_k * betas[int_name]
      base_cohen_d[k] <- MD_k / Sp
    }
    
    for (i in 1:n_boot_samples) {
      set.seed(i)
      sampled_beta <- MASS::mvrnorm(n = 1, mu = betas, Sigma = cov_mat)
      sampled_beta_null <- sampled_beta
      
      for (k in 1:K) {
        int_name <- int_names_k[k]
        if (impact_model %in% c("full", "level")) sampled_beta_null[paste0("indicator_", k)] <- sampled_beta[paste0("indicator_", k)] - betas[paste0("indicator_", k)]
        if (impact_model %in% c("full", "slope")) sampled_beta_null[int_name] <- sampled_beta[int_name] - betas[int_name]
      }
      
      model_copy$coefficients <- sampled_beta
      Sp_boot <- sqrt((var(predict(model_copy, type="response")[data$is_active]) + var(predict(model_copy, type="response", newdata=new_data)[data$is_active])) / 2)
      
      model_copy$coefficients <- sampled_beta_null
      Sp_boot_null <- sqrt((var(predict(model_copy, type="response")[data$is_active]) + var(predict(model_copy, type="response", newdata=new_data)[data$is_active])) / 2)
      
      for (k in 1:K) {
        b_k <- (intervention_lengths[k] - 1) / 2
        ind_name <- paste0("indicator_", k)
        int_name <- int_names_k[k]
        
        if (impact_model == "full") {
          MD_k <- sampled_beta[ind_name] + b_k * sampled_beta[int_name]
          MD_null_k <- sampled_beta_null[ind_name] + b_k * sampled_beta_null[int_name]
        } else if (impact_model == "level") {
          MD_k <- sampled_beta[ind_name]
          MD_null_k <- sampled_beta_null[ind_name]
        } else if (impact_model == "slope") {
          MD_k <- b_k * sampled_beta[int_name]
          MD_null_k <- b_k * sampled_beta_null[int_name]
        }
        cohen_boot_array[i, 1, k] <- MD_k / Sp_boot
        cohen_boot_array[i, 2, k] <- MD_null_k / Sp_boot_null
      }
    }
    for (k in 1:K) ret_vec[k, ] <- c(base_cohen_d[k], quantile(cohen_boot_array[, 1, k], 0.025), quantile(cohen_boot_array[, 1, k], 0.975), mean(abs(cohen_boot_array[, 2, k]) > abs(base_cohen_d[k])))
  }
    
    for (i in 1:n_boot_samples) {
      set.seed(i)
      sampled_beta <- MASS::mvrnorm(n = 1, mu = betas, Sigma = cov_mat)
      sampled_beta_null <- sampled_beta
      for (k in 1:K) {
        if (impact_model %in% c("full", "level")) sampled_beta_null[paste0("indicator_", k)] <- sampled_beta[paste0("indicator_", k)] - betas[paste0("indicator_", k)]
        if (impact_model %in% c("full", "slope")) sampled_beta_null[paste0("shifted_time_", k)] <- sampled_beta[paste0("shifted_time_", k)] - betas[paste0("shifted_time_", k)]
      }
      
      model_copy$coefficients <- sampled_beta
      Sp_boot <- sqrt((var(predict(model_copy, type="response")[data$is_active]) + var(predict(model_copy, type="response", newdata=new_data)[data$is_active])) / 2)
      
      model_copy$coefficients <- sampled_beta_null
      Sp_boot_null <- sqrt((var(predict(model_copy, type="response")[data$is_active]) + var(predict(model_copy, type="response", newdata=new_data)[data$is_active])) / 2)
      
      for (k in 1:K) {
        b_k <- (intervention_lengths[k] - 1) / 2
        ind_name <- paste0("indicator_", k); shift_name <- paste0("shifted_time_", k)
        if (impact_model == "full") {
          MD_k <- sampled_beta[ind_name] + b_k * sampled_beta[shift_name]
          MD_null_k <- sampled_beta_null[ind_name] + b_k * sampled_beta_null[shift_name]
        } else if (impact_model == "level") {
          MD_k <- sampled_beta[ind_name]
          MD_null_k <- sampled_beta_null[ind_name]
        } else if (impact_model == "slope") {
          MD_k <- b_k * sampled_beta[shift_name]
          MD_null_k <- b_k * sampled_beta_null[shift_name]
        }
        cohen_boot_array[i, 1, k] <- MD_k / Sp_boot
        cohen_boot_array[i, 2, k] <- MD_null_k / Sp_boot_null
      }
    }
    for (k in 1:K) ret_vec[k, ] <- c(base_cohen_d[k], quantile(cohen_boot_array[, 1, k], 0.025), quantile(cohen_boot_array[, 1, k], 0.975), mean(abs(cohen_boot_array[, 2, k]) > abs(base_cohen_d[k])))
  }
  
  summar <- summary(model)
  if (isTRUE(print_summary)) print(summar)
  print(ret_vec)
  summar[["Cohen's d"]] <- ret_vec
  if (!isTRUE(counterfactual)) data$predC <- NULL
  
  return(list(model = model, model_summary = summar, data = tibble::as_tibble(data)))
}


#' ITS analysis for continuous outcomes with multiple interruptions
#'
#' \code{its_lm} fits a linear regression model to an Interrupted Time Series with one or more disruptions. 
#' It acts as a unified wrapper, automatically routing the data to either a standard linear model or a 
#' Fourier-adjusted seasonal model based on user input. It returns the fitted model, the summary 
#' (including Cohen's d), and the original data appended with model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the continuous outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right. Should not contain an offset term.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent. Default is NULL.
#' @param multi_model_type A string specifying the multi-disruption architecture. Options are "singular", "multiple", or "mixed". Default is "singular". Note: "mixed" is not supported if seasonality is set to "full".
#' @param freq A positive integer describing the frequency of the time series (e.g., 12 for monthly). Required if seasonality is "full". Default is NULL.
#' @param seasonality A string specifying whether seasonality should be considered. Options are "none" (no seasonal adjustment) or "full" (Fourier terms applied). Default value is "none".
#' @param keep_significant_fourier Logical - if seasonality is "full", indicates whether only significant Fourier terms should be retained. Default is TRUE.
#' @param impact_model A string specifying the assumed impact model: "full", "level", or "slope". Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned. Default is FALSE.
#' @param print_summary Logical - indicating whether the entire model summary should be printed. Default is FALSE.
#' @return The function returns a list with three elements: the fitted linear regression model, the summary of the model (including Cohen's d), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' fit <- its_lm(data=data, form=form, time_name="time", 
#'               intervention_start_indices=start_indices, freq=12, seasonality="none", 
#'               multi_model_type="multiple", impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export

its_lm <- function(data, form, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", freq=NULL, seasonality="none", keep_significant_fourier=TRUE, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  if (!(seasonality %in% c("none", "full"))) {
    stop("Please enter a valid seasonality name. Possible seasonality adjustments are \"none\", and \"full\".")
  }
  
  if (seasonality == "none") {
    mod <- its_lm_wo_seas(data=data, form=form, time_name=time_name, intervention_start_indices=intervention_start_indices, intervention_lengths=intervention_lengths, multi_model_type=multi_model_type, impact_model=impact_model, counterfactual=counterfactual, print_summary=print_summary)
  } else if (seasonality == "full") {
    if (is.null(freq)) {
      stop("Please supply a freq parameter for the seasonality adjustment")
    }
    mod <- its_lm_fourier(data=data, form=form, time_name=time_name, intervention_start_indices=intervention_start_indices, intervention_lengths=intervention_lengths, multi_model_type=multi_model_type, freq=freq, keep_significant_fourier=keep_significant_fourier, impact_model=impact_model, counterfactual=counterfactual, print_summary=print_summary)
  }
  
  return(mod)
}

#' ITS analysis for count outcomes with multiple interruptions
#'
#' \code{its_poisson} acts as a wrapper to fit either a standard or Fourier-adjusted 
#' Poisson (or quasi-Poisson) regression model to an ITS with one or more disruptions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right. Should not contain an offset term.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL. Default value is NULL.
#' @param time_name A string giving the name of the time variable.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent.
#' @param multi_model_type A string specifying the multi-disruption architecture: "singular", "multiple", or "mixed". Default is "singular".
#' @param freq A positive integer describing the frequency of the time series. Required if seasonality is "full".
#' @param seasonality A string specifying whether seasonality should be considered: "none" or "full". Default is "none".
#' @param keep_significant_fourier Logical - if seasonality is "full", indicates whether only significant Fourier terms should be retained. Default is TRUE.
#' @param over_dispersion Logical - if TRUE, uses quasi-Poisson. Default is FALSE. (Note: ignored for mixed-effects models).
#' @param impact_model A string specifying the impact model: "full", "level", or "slope". Default is "full".
#' @param counterfactual Logical - indicates whether model-based counterfactual values should be returned. Default is FALSE.
#' @param print_summary Logical - indicates whether the model summary should be printed. Default is FALSE.
#' @return A list with three elements: the fitted model, the model summary (including relative risk), and the data with predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' fit <- its_poisson(data=data, form=form, offset_name="labour", time_name="time", 
#'                    intervention_start_indices=start_indices, multi_model_type="multiple", 
#'                    seasonality="none", impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export

its_poisson <- function(data, form, offset_name=NULL, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", freq=NULL, seasonality="none", keep_significant_fourier=TRUE, over_dispersion=FALSE, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  if (!(seasonality %in% c("none", "full"))) {
    stop("Please enter a valid seasonality name. Possible seasonality adjustments are \"none\", and \"full\".")
  }
  
  if (seasonality == "none") {
    mod <- its_poisson_wo_seas(data=data, form=form, offset_name=offset_name, time_name=time_name, intervention_start_indices=intervention_start_indices, intervention_lengths=intervention_lengths, multi_model_type=multi_model_type, over_dispersion=over_dispersion, impact_model=impact_model, counterfactual=counterfactual, print_summary=print_summary)
  } else if (seasonality == "full") {
    if (is.null(freq)) {
      stop("Please supply a freq parameter for the seasonality adjustment")
    }
    mod <- its_poisson_fourier(data=data, form=form, offset_name=offset_name, time_name=time_name, intervention_start_indices=intervention_start_indices, intervention_lengths=intervention_lengths, multi_model_type=multi_model_type, over_dispersion=over_dispersion, freq=freq, keep_significant_fourier=keep_significant_fourier, impact_model=impact_model, counterfactual=counterfactual, print_summary=print_summary)
  }
  
  return(mod)
}

#' ITS analysis for zero-inflated count outcomes with multiple interruptions
#'
#' \code{its_zero_inflated_poisson} acts as a wrapper to fit a zero-inflated Poisson 
#' regression model to an ITS with one or more disruptions, with optional seasonal adjustment.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables.
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL.
#' @param time_name A string giving the name of the time variable.
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention.
#' @param multi_model_type A string specifying the multi-disruption architecture: "singular", "multiple", or "mixed". Default is "singular".
#' @param freq A positive integer describing the frequency of the time series. Required if seasonality is TRUE.
#' @param seasonality Logical - indicating whether seasonal adjustment via Fourier terms should be used. Default is FALSE.
#' @param impact_model A string specifying the impact model: "full", "level", or "slope". Default is "full".
#' @param counterfactual Logical - indicating whether the counterfactual values should be returned. Default is FALSE.
#' @param print_summary Logical - indicating whether the model summary should be printed. Default is FALSE.
#' @return A list with three elements: the fitted zero-inflated model, the model summary (including relative risk), and the data with predictions.
#' @examples
#' data <- zero_inflated_sim_data
#' form <- as.formula("monthly_total ~ time")
#' start_indices <- c(which(data$Year==2020 & data$Month==3), which(data$Year==2021 & data$Month==1))
#' fit <- its_zero_inflated_poisson(data=data, form=form, time_name="time", 
#'                                  intervention_start_indices=start_indices, freq=12, 
#'                                  seasonality=TRUE, impact_model="full", counterfactual=TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @importFrom pscl zeroinfl
#' @export

its_zero_inflated_poisson <- function(data, form, offset_name=NULL, time_name, intervention_start_indices, intervention_lengths=NULL, multi_model_type="singular", freq=NULL, seasonality=FALSE, impact_model="full", counterfactual=FALSE, print_summary=FALSE) {
  
  if (!isTRUE(seasonality)) {
    out <- its_zero_inflated(data = data, form = form, offset_name = offset_name, time_name = time_name, intervention_start_indices = intervention_start_indices, intervention_lengths = intervention_lengths, multi_model_type = multi_model_type, impact_model = impact_model, counterfactual = counterfactual, print_summary = print_summary)
  } else {
    if (is.null(freq)) {
      stop("Please supply a freq parameter for the seasonality adjustment")
    }
    out <- its_zero_inflated_fourier(data = data, form = form, offset_name = offset_name, time_name = time_name, intervention_start_indices = intervention_start_indices, intervention_lengths = intervention_lengths, multi_model_type = multi_model_type, freq = freq, impact_model = impact_model, counterfactual = counterfactual, print_summary = print_summary)
  }
  
  return(out)
}


#' Plot ITS fitted values with multiple interruptions
#'
#' \code{plot_its_lm} uses ggplot2 to plot the model-based fitted values, together with a scatterplot 
#' of the observed time series. It supports multiple disruptions, automatically drawing broken regression 
#' lines and vertical segments at the start and end of each interruption.
#'
#' @param data The data frame containing the original data and the model predictions (returned by the its_lm function).
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent.
#' @param y_lab A string with the y-axis label for the plot.
#' @param response A string giving the name of the response variable to be plotted in the scatterplot.
#' @param date_name A string giving the name of the date column. The date column must be a Date object.
#' @return A ggplot object including a scatterplot of the time series, the predictions line, and the counterfactual predictions in red (if available).
#' @examples
#' \dontrun{
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' fit <- its_lm(data=data, form=form, time_name="time", 
#'               intervention_start_indices=start_indices, 
#'               multi_model_type="multiple", counterfactual=TRUE)
#' 
#' p <- plot_its_lm(data=fit$data, intervention_start_indices=start_indices, 
#'                  y_lab="Unemployment percent", response="percent", date_name="dt")
#' print(p)
#' }
#' @importFrom lubridate is.Date
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_point ylab xlab geom_line draw_key_smooth geom_segment scale_color_manual scale_linetype_manual theme element_blank theme_bw
#' @export

plot_its_lm <- function(data, intervention_start_indices, intervention_lengths = NULL, y_lab, response, date_name){
  pred <- predC <- indicator <- NULL
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(response %in% colnames(data))) stop("Please make sure that the response variable belongs to colnames(data)")
  if (!(date_name %in% colnames(data))) stop("Please make sure that the date_name column belongs to colnames(data)")
  else if (!(lubridate::is.Date(data[[date_name]]))) stop("Please make sure that the date_name column is a Date object")
  
  n <- nrow(data)
  if (!all(intervention_start_indices > 1 & intervention_start_indices <= n)) stop("Please make sure that all intervention_start_indices are valid.")
  if (!("pred" %in% colnames(data))) stop("The data needs to include the column of predictions.")
  if (!(is.character(y_lab))) { y_lab <- ""; print("y_lab is not a string and thus will not by used") }
  
  # BOUNDARY CALCULATION 
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  boundaries <- intervention_start_indices
  for (k in 1:length(intervention_start_indices)) {
    end_idx <- intervention_start_indices[k] + intervention_lengths[k]
    if (end_idx <= n) boundaries <- c(boundaries, end_idx)
  }
  boundaries <- sort(unique(boundaries))
  boundaries <- boundaries[boundaries > 1 & boundaries <= n]
  
  data$group_id <- 1
  for (idx in boundaries) data$group_id[idx:n] <- data$group_id[idx:n] + 1
  
  # Safe string extraction for boundaries
  segments_df <- data.frame(
    xdot = data[[date_name]][boundaries],
    y_start = data$pred[boundaries - 1],
    y_end = data$pred[boundaries]
  )
  
  # Convert to symbols ONLY for ggplot
  resp_sym <- rlang::sym(response)
  date_sym <- rlang::sym(date_name)
  
  # PLOTTING 
  if ("predC" %in% colnames(data)){
    Group <- levels( factor(c("Counterfactual","Fitted values")))
    p <- ggplot2::ggplot(data = data , ggplot2::aes(y = !!resp_sym, x = !!date_sym)) +
      ggplot2::geom_point()  + ggplot2::ylab(y_lab) + ggplot2::xlab("Date") +
      ggplot2::geom_line(ggplot2::aes(y = predC, x = !!date_sym, color="Counterfactual", linetype = "Counterfactual"), linewidth=1, key_glyph = ggplot2::draw_key_smooth, na.rm=TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = pred, x = !!date_sym, group = group_id, color="Fitted values", linetype = "Fitted values"), linewidth=1, key_glyph = ggplot2::draw_key_smooth) +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x=xdot, xend=xdot, y=y_start, yend=y_end), inherit.aes = FALSE, linetype=3, linewidth=1) +
      ggplot2::scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"), breaks = Group) +
      ggplot2::scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1), breaks = Group) +
      ggplot2::theme_bw() + ggplot2::theme(legend.key = ggplot2::element_blank())
  } else {
    p <- ggplot2::ggplot(data = data , ggplot2::aes(y = !!resp_sym, x = !!date_sym)) +
      ggplot2::geom_point()  + ggplot2::ylab(y_lab) + ggplot2::xlab("Date") +
      ggplot2::geom_line(ggplot2::aes(y = pred, x = !!date_sym, group = group_id), linewidth=1, key_glyph = ggplot2::draw_key_smooth) +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x=xdot, xend=xdot, y=y_start, yend=y_end), inherit.aes = FALSE, linetype=3, linewidth=1) +
      ggplot2::theme_bw() + ggplot2::theme(legend.key = ggplot2::element_blank())
  }
  return(p)
}


plot_its_poisson <- function(data, intervention_start_indices, intervention_lengths = NULL, y_lab, response, offset_name = NULL, date_name){
  pred <- predC <- indicator <- NULL
  if (!(is.data.frame(data) | tibble::is_tibble(data))) stop("Please make sure that data is either a data frame or a tibble")
  if (!(response %in% colnames(data))) stop("Please make sure that the response variable belongs to colnames(data)")
  if (!(date_name %in% colnames(data))) stop("Please make sure that the date_name column belongs to colnames(data)")
  else if (!(lubridate::is.Date(data[[date_name]]))) stop("Please make sure that the date_name column is a Date object")
  
  n <- nrow(data)
  if (!all(intervention_start_indices > 1 & intervention_start_indices <= n)) stop("Please make sure all intervention_start_indices are valid.")
  if (!("pred" %in% colnames(data))) stop("The data needs to include the column of predictions.")
  if (!(is.character(y_lab))) { y_lab <- ""; print("y_lab is not a string and thus will not by used") }
  
  if (!(is.null(offset_name))){
    if (!(offset_name %in% colnames(data))) stop("Please make sure that offset_name belongs to colnames(data)")
    else {
      data$pred <- data$pred * 100 / data[[offset_name]]
      data[[response]] <- data[[response]] * 100 / data[[offset_name]]
      if ("predC" %in% colnames(data)) data$predC <- data$predC * 100 / data[[offset_name]]
    }
  }
  
  # BOUNDARY CALCULATION 
  if (is.null(intervention_lengths)) intervention_lengths <- n - intervention_start_indices + 1
  boundaries <- intervention_start_indices
  for (k in 1:length(intervention_start_indices)) {
    end_idx <- intervention_start_indices[k] + intervention_lengths[k]
    if (end_idx <= n) boundaries <- c(boundaries, end_idx)
  }
  boundaries <- sort(unique(boundaries))
  boundaries <- boundaries[boundaries > 1 & boundaries <= n]
  
  data$group_id <- 1
  for (idx in boundaries) data$group_id[idx:n] <- data$group_id[idx:n] + 1
  
  # Safe string extraction for boundaries
  segments_df <- data.frame(
    xdot = data[[date_name]][boundaries],
    y_start = data$pred[boundaries - 1],
    y_end = data$pred[boundaries]
  )
  
  # Convert to symbols ONLY for ggplot
  resp_sym <- rlang::sym(response)
  date_sym <- rlang::sym(date_name)
  
  # PLOTTING 
  if ("predC" %in% colnames(data)){
    Group <- levels( factor(c("Counterfactual","Fitted values")))
    p <- ggplot2::ggplot(data = data , ggplot2::aes(y = !!resp_sym, x = !!date_sym)) +
      ggplot2::geom_point()  + ggplot2::ylab(y_lab) + ggplot2::xlab("Date") +
      ggplot2::geom_line(ggplot2::aes(y = predC, x = !!date_sym, color="Counterfactual", linetype = "Counterfactual"), linewidth=1, key_glyph = ggplot2::draw_key_smooth, na.rm=TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = pred, x = !!date_sym, group = group_id, color="Fitted values", linetype = "Fitted values"), linewidth=1, key_glyph = ggplot2::draw_key_smooth) +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x=xdot, xend=xdot, y=y_start, yend=y_end), inherit.aes = FALSE, linetype=3, linewidth=1) +
      ggplot2::scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"), breaks = Group) +
      ggplot2::scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1), breaks = Group) +
      ggplot2::theme_bw() + ggplot2::theme(legend.key = ggplot2::element_blank())
  } else {
    p <- ggplot2::ggplot(data = data , ggplot2::aes(y = !!resp_sym, x = !!date_sym)) +
      ggplot2::geom_point()  + ggplot2::ylab(y_lab) + ggplot2::xlab("Date") +
      ggplot2::geom_line(ggplot2::aes(y = pred, x = !!date_sym, group = group_id), linewidth=1, key_glyph = ggplot2::draw_key_smooth) +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x=xdot, xend=xdot, y=y_start, yend=y_end), inherit.aes = FALSE, linetype=3, linewidth=1) +
      ggplot2::theme_bw() + ggplot2::theme(legend.key = ggplot2::element_blank())
  }
  return(p)
}

#' Plot ITS fitted values for count outcomes with multiple interruptions
#'
#' \code{plot_its_poisson} uses ggplot2 to plot the model-based fitted values, together with a scatterplot 
#' of the observed time series. If an offset is provided, it automatically normalizes the predictions 
#' and observations into a standardized rate. It supports multiple disruptions by drawing broken 
#' regression lines and vertical segments at the boundaries of each interruption.
#'
#' @param data The data frame containing the original data and the model predictions (returned by its_poisson or its_zero_inflated_poisson).
#' @param intervention_start_indices Numeric vector - stating the time point indices of the start of each intervention.
#' @param intervention_lengths Numeric vector - stating the duration of each intervention. If NULL, disruptions are assumed to be permanent.
#' @param y_lab A string with the y-axis label for the plot.
#' @param response A string giving the name of the response variable to be plotted in the scatterplot.
#' @param offset_name A string indicating the name of the offset column in the data, or NULL. If provided, values are normalized as (value * 100 / offset).
#' @param date_name A string giving the name of the date column. The date column must be a Date object.
#' @return A ggplot object including a scatterplot of the time series, the predictions line, and the counterfactual predictions in red (if available).
#' @examples
#' \dontrun{
#' data <- unemployed 
#' form <- as.formula("unemployed ~ time")
#' start_indices <- c(which(data$year==2020 & data$month==3), which(data$year==2021 & data$month==1))
#' fit <- its_poisson(data=data, form=form, offset_name="labour", time_name="time", 
#'                    intervention_start_indices=start_indices, 
#'                    multi_model_type="multiple", counterfactual=TRUE)
#' 
#' p <- plot_its_poisson(data=fit$data, intervention_start_indices=start_indices, 
#'                       y_lab="Unemployment Rate", response="unemployed", 
#'                       offset_name="labour", date_name="dt")
#' print(p)
#' }
#' @importFrom lubridate is.Date
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_point ylab xlab geom_line draw_key_smooth geom_segment scale_color_manual scale_linetype_manual theme element_blank theme_bw
#' @export

plot_its_poisson <- function(data, intervention_start_indices, intervention_lengths = NULL, y_lab, response, offset_name = NULL, date_name){
  pred <- predC <- indicator <- NULL
  if (!(is.data.frame(data) | tibble::is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable belongs to colnames(data)")
  }
  if (!(date_name %in% colnames(data))){
    stop("Please make sure that the date_name column belongs to colnames(data)")
  } else if (!(lubridate::is.Date(data[[date_name]]))){
    stop("Please make sure that the date_name column is a Date object")
  }
  
  n <- nrow(data)
  if (!all(intervention_start_indices > 1 & intervention_start_indices <= n)){
    stop("Please make sure that all intervention_start_indices are greater than 1 and less or equal to nrow(data).")
  }
  if (!("pred" %in% colnames(data))){
    stop("The data needs to include the column of predictions. Please use as input the data output of the function its_poisson()")
  }
  if (!(is.character(y_lab))){
    y_lab <- ""
    print("y_lab is not a string and thus will not by used")
  }
  
  # Handle Offset before creating plotting boundaries
  if (!(is.null(offset_name))){
    if (!(offset_name %in% colnames(data))){
      stop("Please make sure that offset_name belongs to colnames(data)")
    } else {
      data$pred <- data$pred * 100 / data[[offset_name]]
      data[[response]] <- data[[response]] * 100 / data[[offset_name]]
      if ("predC" %in% colnames(data)){
        data$predC <- data$predC * 100 / data[[offset_name]]
      }
    }
  }
  
  # BOUNDARY CALCULATION 
  if (is.null(intervention_lengths)) {
    intervention_lengths <- n - intervention_start_indices + 1
  }
  
  boundaries <- intervention_start_indices
  for (k in 1:length(intervention_start_indices)) {
    end_idx <- intervention_start_indices[k] + intervention_lengths[k]
    if (end_idx <= n) {
      boundaries <- c(boundaries, end_idx)
    }
  }
  boundaries <- sort(unique(boundaries))
  boundaries <- boundaries[boundaries > 1 & boundaries <= n]
  
  # Create a grouping ID to elegantly break the geom_line at every boundary
  data$group_id <- 1
  for (idx in boundaries) {
    data$group_id[idx:n] <- data$group_id[idx:n] + 1
  }
  
  # Create a dataframe for the vertical dashed segments for ALL boundaries
  segments_df <- data.frame(
    xdot = data[[date_name]][boundaries],
    y_start = data$pred[boundaries - 1],
    y_end = data$pred[boundaries]
  )
  
  response <- rlang::sym(response)
  date_name <- rlang::sym(date_name)
  
  # PLOTTING 
  if ("predC" %in% colnames(data)){
    Group <- levels( factor(c("Counterfactual","Fitted values")))
    p <- ggplot2::ggplot(data = data , ggplot2::aes(y = !!response, x = !!date_name)) +
      ggplot2::geom_point()  + ggplot2::ylab(y_lab) + ggplot2::xlab("Date") +
      ggplot2::geom_line(ggplot2::aes(y = predC, x = !!date_name, color="Counterfactual", linetype = "Counterfactual"), linewidth=1, key_glyph = ggplot2::draw_key_smooth, na.rm=TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = pred, x = !!date_name, group = group_id, color="Fitted values", linetype = "Fitted values"), linewidth=1, key_glyph = ggplot2::draw_key_smooth) +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x=xdot, xend=xdot, y=y_start, yend=y_end), inherit.aes = FALSE, linetype=3, linewidth=1) +
      ggplot2::scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"), breaks = Group) +
      ggplot2::scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1), breaks = Group) +
      ggplot2::theme_bw() + ggplot2::theme(legend.key = ggplot2::element_blank())
  } else {
    p <- ggplot2::ggplot(data = data , ggplot2::aes(y = !!response, x = !!date_name)) +
      ggplot2::geom_point()  + ggplot2::ylab(y_lab) + ggplot2::xlab("Date") +
      ggplot2::geom_line(ggplot2::aes(y = pred, x = !!date_name, group = group_id), linewidth=1, key_glyph = ggplot2::draw_key_smooth) +
      ggplot2::geom_segment(data = segments_df, ggplot2::aes(x=xdot, xend=xdot, y=y_start, yend=y_end), inherit.aes = FALSE, linetype=3, linewidth=1) +
      ggplot2::theme_bw() + ggplot2::theme(legend.key = ggplot2::element_blank())
  }
  
  return(p)
}
