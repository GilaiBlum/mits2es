#' Build Disruption Features for Multiple ITS Models
#'
#' This helper function takes a dataframe and user inputs for multiple disruptions,
#' returning the dataframe appended with K level indicators and K shifted time continuous vectors.
#'
#' @param data A data frame or tibble containing the time series.
#' @param time_name A string giving the name of the time variable in `data`.
#' @param start_indices A numeric vector of length K stating the start time indices of each intervention.
#' @param lengths A numeric vector of length K stating the duration of each intervention. 
#'                If NULL, disruptions are assumed to be permanent (lasting until the end of the series).
#' @return The original data frame with 2*K new columns appended (indicator_k and shifted_time_k).
build_disruption_features <- function(data, time_name, start_indices, lengths = NULL) {
  
  # Basic validation
  if (!is.data.frame(data) && !tibble::is_tibble(data)) {
    stop("Please make sure that data is either a data frame or a tibble.")
  }
  if (!(time_name %in% colnames(data))) {
    stop("Please make sure that time_name belongs to colnames(data).")
  }
  
  n <- nrow(data)
  K <- length(start_indices)
  
  # Handle disruption lengths (defaults to permanent if NULL)
  if (is.null(lengths)) {
    lengths <- n - start_indices + 1
  } else if (length(lengths) != K) {
    stop("The length of 'lengths' must match the length of 'start_indices'.")
  }
  
  # Generate matrices for each K
  for (k in 1:K) {
    start_ind <- start_indices[k]
    len <- lengths[k]
    
    # Ensure temporary disruptions do not exceed the dataframe's boundaries
    end_ind <- min(n, start_ind + len - 1) 
    
    if (start_ind > 1 && start_ind <= n) {
      
      # 1. Level Indicator Matrix (I_k)
      ind_col <- paste0("indicator_", k)
      data[[ind_col]] <- 0
      data[[ind_col]][start_ind:end_ind] <- 1
      
      # 2. Shifted Time Matrix (t - t_k^*)
      # Following the package's original architecture, this calculates the time shift globally.
      # The interaction term in the formula (indicator_k:shifted_time_k) will zero it out pre-intervention.
      shift_col <- paste0("shifted_time_", k)
      data[[shift_col]] <- data[[time_name]] - start_ind
      
    } else {
      warning(paste("Disruption", k, "start index is out of bounds. Skipping."))
    }
  }
  
  return(data)
}

#' Build Disruption Features for Mixed-Effects ITS
#'
#' This helper generates unified feature columns and a grouping factor (disruption_id)
#' required by lme4, ensuring the time series degrees of freedom are preserved.
#'
#' @param data A data frame or tibble.
#' @param time_name A string giving the name of the time variable.
#' @param start_indices A numeric vector of length K stating the start time indices.
#' @param lengths A numeric vector of length K stating durations.
#' @return The original data frame with 3 new columns: generic_indicator, generic_shifted_time, disruption_id.
build_mixed_features <- function(data, time_name, start_indices, lengths = NULL) {
  
  n <- nrow(data)
  K <- length(start_indices)
  
  if (is.null(lengths)) {
    lengths <- n - start_indices + 1
  }
  
  # Initialize base state
  data$generic_indicator <- 0
  data$generic_shifted_time <- 0
  data$disruption_id <- "Disruption_1"
  
  for (k in 1:K) {
    start_ind <- start_indices[k]
    len <- lengths[k]
    end_ind <- min(n, start_ind + len - 1)
    
    if (start_ind > 1 && start_ind <= n) {
      
      # Safety Check for Overlapping Disruptions
      if (any(data$generic_indicator[start_ind:end_ind] == 1)) {
        warning(paste("Disruption", k, "overlaps with a previous disruption. The unified Mixed-Effects architecture assumes sequential events. Overlapping days will be reassigned to Disruption", k))
      }
      
      data$generic_indicator[start_ind:end_ind] <- 1
      data$generic_shifted_time[start_ind:end_ind] <- data[[time_name]][start_ind:end_ind] - start_ind
      data$disruption_id[start_ind:end_ind] <- paste0("Disruption_", k)
    }
  }
  
  # Convert to factor for lme4 grouping
  data$disruption_id <- as.factor(data$disruption_id)
  
  return(data)
}