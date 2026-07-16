#!/usr/bin/env Rscript
# Run one structured robustness scenario and write a grid-compatible RDS
# (class_metrics + support_contrast_summary) for sim_study_robustness.R.
#
# Env (set by sim_study_robustness.R run_one):
#   PP_SCENARIO_FAMILY   effect_modification | geometry_transport
#   PP_SCENARIO_ID
#   PP_OUTPUT_BASENAME
#   PP_H_TRUE            (effect_modification)
#   PP_GEOMETRY_M        (geometry_transport)
#   PP_SIMS, PP_TARGET_POINTS, PP_CONTROL_K, PP_TREATED_K
#   PP_SEM_INNER_ITER, PP_SEM_OUTER_ITER (optional)
#   PP_TRUTH_SIMS, PP_FORWARD_SIMS (optional)

suppressPackageStartupMessages({
  library(PPDisentangle)
  library(spatstat)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) < 1L || i[[1L]] >= length(args)) return(default)
  args[[i[[1L]] + 1L]]
}
arg_int <- function(flag, default, minimum = 1L) {
  x <- suppressWarnings(as.integer(get_arg(flag, default)))
  if (!is.finite(x) || is.na(x)) x <- default
  max(as.integer(minimum), x)
}
env_num <- function(name, default) {
  x <- suppressWarnings(as.numeric(Sys.getenv(name, "")))
  if (!is.finite(x) || is.na(x)) default else x
}
env_int <- function(name, default, minimum = 1L) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, "")))
  if (!is.finite(x) || is.na(x)) x <- default
  max(as.integer(minimum), x)
}

script_dir <- {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1L]])))
  else normalizePath(file.path(getwd(), "inst", "sim_study"), mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(script_dir, "structured_robustness_utils.R"), local = FALSE)

test_mode <- "--test" %in% args
family <- Sys.getenv("PP_SCENARIO_FAMILY", "")
scenario_id <- Sys.getenv("PP_SCENARIO_ID", family)
run_basename <- Sys.getenv("PP_OUTPUT_BASENAME", paste0("structured_", family))
if (!nzchar(family) || !family %in% c("effect_modification", "geometry_transport")) {
  stop("PP_SCENARIO_FAMILY must be effect_modification or geometry_transport")
}

n_rep <- arg_int("--sims", env_int("PP_SIMS", if (test_mode) 1L else 32L))
if (test_mode) n_rep <- min(n_rep, 1L)
n_truth <- env_int(
  "PP_TRUTH_SIMS",
  env_int("PP_STRUCTURED_TRUTH_SIMS", if (test_mode) 2L else 500L)
)
n_cf <- env_int(
  "PP_FORWARD_SIMS",
  env_int("PP_STRUCTURED_FORWARD_SIMS", if (test_mode) 2L else 100L)
)
sem_inner <- env_int("PP_SEM_INNER_ITER", if (test_mode) 5L else 1000L)
sem_outer <- env_int("PP_SEM_OUTER_ITER", if (test_mode) 1L else 3L)
sem_labels <- env_int("PP_SEM_N_LABELLINGS", if (test_mode) 1L else 10L, 0L)
sem_props <- env_int("PP_SEM_N_PROPS", if (test_mode) 1L else 20L)
fit_maxit <- if (test_mode) 30L else 1000L
base_seed <- env_int("PP_STRUCTURED_SEED", 20260714L)
h_true <- env_num("PP_H_TRUE", 0)
geometry_m <- suppressWarnings(as.integer(Sys.getenv("PP_GEOMETRY_M", "NA")))
control_K <- env_num("PP_CONTROL_K", 0.8)
treated_K <- env_num("PP_TREATED_K", 0.2)
target_points <- env_num("PP_TARGET_POINTS", 2500)
rep_workers <- min(
  n_rep,
  env_int("PP_STRUCTURED_REP_WORKERS", env_int("PP_CPUS", 1L))
)
if (test_mode || identical(.Platform$OS.type, "windows")) rep_workers <- 1L

structured_lapply <- function(X, FUN) {
  if (rep_workers <= 1L || length(X) <= 1L) return(lapply(X, FUN))
  parallel::mclapply(
    X, FUN,
    mc.cores = min(rep_workers, length(X)),
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

omega <- spatstat.geom::owin(c(0, 100), c(0, 100))
partition <- spatstat.geom::quadrats(X = omega, nx = 10, ny = 10)
treatment_time <- 10
end_time <- 110
post_window <- c(treatment_time, end_time)
expected_per_mu <- treatment_time / (1 - control_K) +
  (end_time - treatment_time) * (0.5 / (1 - control_K) + 0.5 / (1 - treated_K))
true_mu <- target_points / expected_per_mu
control_params <- list(
  mu = true_mu, alpha = 0.01, beta = 10, K = control_K,
  kernel = "exponential", spatial_kernel = "exponential"
)
treated_params <- list(
  mu = true_mu, alpha = 0.01, beta = 10, K = treated_K,
  kernel = "exponential", spatial_kernel = "exponential"
)
sem_control <- list(
  N_labellings = sem_labels,
  N_iter = sem_outer,
  iter = sem_inner,
  n_props = sem_props,
  change_factor = 0.005,
  param_update_cadence = env_int("PP_SEM_PARAM_UPDATE_CADENCE", 10L),
  proposal_update_cadence = env_int("PP_SEM_PROPOSAL_UPDATE_CADENCE", 1L),
  update_control_params = FALSE,
  include_starting_data = TRUE,
  update_starting_data = TRUE,
  maxit = fit_maxit,
  verbose = FALSE
)

simulate_pre <- function(rep_id) {
  set.seed(structured_stage_seed(base_seed, 10L, rep_id))
  pre <- PPDisentangle::sim_hawkes(
    params = control_params,
    windowT = c(0, treatment_time),
    windowS = omega,
    optimized = TRUE
  )
  prepare_pre_history(data.frame(
    x = pre$x, y = pre$y, t = pre$t, n = length(pre$t),
    background = pre$background
  ), partition)
}

bind_catalogue_local <- function(pre, post, grid, z, h) {
  pre <- attach_design_columns(pre, partition, grid, z, h = 0)
  post <- attach_design_columns(post, partition, grid, z, h = h)
  # Drop pre-only inferred_process so bind_rows does not create NA rows that
  # break adaptive_SEM (it only auto-fills when the column is absent).
  pre$inferred_process <- NULL
  cat <- dplyr::bind_rows(pre, post)
  if (!"inferred_process" %in% names(cat) || anyNA(cat$inferred_process)) {
    cat$inferred_process <- ifelse(
      cat$t < treatment_time, "control", as.character(cat$location_process)
    )
  }
  cat
}

# Local copy: structured script defines this helper; keep scenario runner self-contained.
contrast_effect_regimes <- function(regime_df) {
  need <- c("z_flip_plus", "z_flip_minus", "z_all", "z_none")
  if (!all(need %in% regime_df$allocation)) {
    stop(
      "contrast_effect_regimes expects allocations: ",
      paste(need, collapse = ", ")
    )
  }
  getm <- function(a) regime_df$mean[regime_df$allocation == a]
  getse <- function(a) regime_df$mc_se[regime_df$allocation == a]
  data.frame(
    target = c("psi_global", "psi_flip_plus", "psi_flip_minus"),
    estimate = c(
      getm("z_all") - getm("z_none"),
      getm("z_flip_plus") - getm("z_none"),
      getm("z_flip_minus") - getm("z_none")
    ),
    mc_se = c(
      sqrt(getse("z_all")^2 + getse("z_none")^2),
      sqrt(getse("z_flip_plus")^2 + getse("z_none")^2),
      sqrt(getse("z_flip_minus")^2 + getse("z_none")^2)
    ),
    stringsAsFactors = FALSE
  )
}

label_row_to_recovery <- function(lr_all, method, sim_id) {
  d <- lr_all[lr_all$stratum == "all", , drop = FALSE]
  if (nrow(d) < 1L) return(NULL)
  data.frame(
    method = method,
    sim_id = as.integer(sim_id),
    n_post = as.integer(d$n[[1L]]),
    accuracy = as.numeric(d$accuracy[[1L]]),
    balanced_accuracy = as.numeric(d$balanced_accuracy[[1L]]),
    precision_treated = NA_real_,
    recall_treated = as.numeric(d$recall_treated[[1L]]),
    specificity_control = as.numeric(d$recall_control[[1L]]),
    f1_treated = NA_real_,
    stringsAsFactors = FALSE
  )
}

summarize_class <- function(lr_df) {
  if (is.null(lr_df) || nrow(lr_df) < 1L) return(NULL)
  lr_df %>%
    mutate(method = ifelse(.data$method == "SEM_adaptive", "SEM", .data$method)) %>%
    group_by(.data$method) %>%
    summarize(
      n_kept = n(),
      mean_accuracy = mean(.data$accuracy, na.rm = TRUE),
      mean_balanced_accuracy = mean(.data$balanced_accuracy, na.rm = TRUE),
      mean_precision_treated = mean(.data$precision_treated, na.rm = TRUE),
      mean_recall_treated = mean(.data$recall_treated, na.rm = TRUE),
      mean_specificity_control = mean(.data$specificity_control, na.rm = TRUE),
      mean_f1_treated = mean(.data$f1_treated, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_support <- function(psi_df) {
  if (is.null(psi_df) || nrow(psi_df) < 1L) return(NULL)
  psi_df %>%
    mutate(
      labelling = ifelse(.data$method == "SEM_adaptive", "SEM", .data$method),
      abs_error = abs(.data$psi_estimate - .data$psi_truth),
      pct_error = ifelse(
        is.finite(.data$psi_truth) & abs(.data$psi_truth) > 1e-8,
        100 * .data$abs_error / abs(.data$psi_truth),
        NA_real_
      ),
      bias = .data$psi_estimate - .data$psi_truth
    ) %>%
    group_by(.data$labelling, .data$contrast_family) %>%
    summarize(
      n_total = n(),
      n = sum(is.finite(.data$psi_estimate) & is.finite(.data$psi_truth)),
      n_failed = sum(!is.finite(.data$psi_estimate) | !is.finite(.data$psi_truth)),
      mean_psi_estimate = mean(.data$psi_estimate, na.rm = TRUE),
      mean_psi_truth = mean(.data$psi_truth, na.rm = TRUE),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      mean_pct_error = mean(.data$pct_error, na.rm = TRUE),
      mean_bias = mean(.data$bias, na.rm = TRUE),
      se_bias = {
        ok <- is.finite(.data$bias)
        if (sum(ok) > 1L) stats::sd(.data$bias[ok]) / sqrt(sum(ok)) else NA_real_
      },
      .groups = "drop"
    )
}

message(sprintf(
  paste0(
    "[structured-scenario] family=%s id=%s sims=%d rep_workers=%d ",
    "sem_inner=%d truth_sims=%d forward_sims=%d"
  ),
  family, scenario_id, n_rep, rep_workers, sem_inner, n_truth, n_cf
))

failed_fit <- function(message, h = NA_real_) {
  list(
    par = c(mu = NA_real_, alpha = NA_real_, beta = NA_real_,
            K = NA_real_, h = h),
    converged = FALSE,
    convergence = 99L,
    message = message
  )
}

fit_status_row <- function(fit, method, rep_id) {
  data.frame(
    method = method,
    sim_id = as.integer(rep_id),
    converged = isTRUE(fit$converged),
    convergence = as.integer(fit$convergence %||% if (isTRUE(fit$converged)) 0L else 99L),
    message = as.character(fit$message %||% ""),
    stringsAsFactors = FALSE
  )
}

collect_replication_parts <- function(parts, family_label) {
  ok <- vapply(parts, function(x) is.list(x) && isTRUE(x$ok), logical(1L))
  if (!all(ok)) {
    failures <- vapply(parts[!ok], function(x) {
      if (is.list(x) && nzchar(x$error %||% "")) x$error else "worker returned no result"
    }, character(1L))
    stop(sprintf(
      "%s failed in %d/%d replications: %s",
      family_label, sum(!ok), length(parts),
      paste(utils::head(failures, 3L), collapse = " | ")
    ))
  }
  list(
    labels = bind_rows(lapply(parts, `[[`, "labels")),
    psi = bind_rows(lapply(parts, `[[`, "psi")),
    fit_status = bind_rows(lapply(parts, `[[`, "fit_status"))
  )
}

assert_complete_rows <- function(label_df, psi_df, fit_status, contrasts) {
  for (method in c("naive", "SEM")) {
    n_label <- sum(label_df$method == method)
    n_fit <- sum(fit_status$method == method)
    if (n_label != n_rep || n_fit != n_rep) {
      stop(sprintf(
        "Incomplete %s outputs: labels=%d fits=%d expected=%d",
        method, n_label, n_fit, n_rep
      ))
    }
    for (contrast in contrasts) {
      n_psi <- sum(psi_df$method == method & psi_df$contrast_family == contrast)
      if (n_psi != n_rep) {
        stop(sprintf(
          "Incomplete %s/%s causal outputs: %d expected=%d",
          method, contrast, n_psi, n_rep
        ))
      }
    }
  }
}

label_recovery_prefilter <- NULL
psi_rows <- NULL
fit_status <- NULL
coarseness_val <- NA_real_

if (identical(family, "effect_modification")) {
  design <- make_effect_modification_design(partition, seed = base_seed + 101L)
  h_level <- match(round(h_true, 8L), round(c(0, 0.3, 0.6), 8L))
  if (is.na(h_level)) stop("Unsupported PP_H_TRUE: ", h_true)
  run_effect_replication <- function(rep_id) {
    tryCatch({
    label_rows <- list()
    causal_rows <- list()
    message(sprintf("[structured-scenario] effect h=%.3f rep %d/%d", h_true, rep_id, n_rep))
    pre <- simulate_pre(rep_id + 1000L * h_level)
    set.seed(structured_stage_seed(base_seed, 20L, rep_id, round(100 * h_true)))
    post <- simulate_structured_catalogue(
      partition, omega, post_window, design$z_obs,
      control_params, treated_params, pre, design$grid, h_true
    )
    catalogue <- bind_catalogue_local(pre, post, design$grid, design$z_obs, h_true)
    labels <- make_oracle_and_naive_labels(catalogue, treatment_time)
    sem <- run_effect_aware_sem(
      catalogue, partition, design$grid, design$z_obs, omega,
      post_window, treatment_time, control_params, treated_params,
      h_init = 0, interaction_rounds = 1L,
      sem_control = sem_control,
      seed = structured_stage_seed(base_seed, 21L, rep_id, round(100 * h_true))
    )
    if (is.null(sem$labelling)) stop("Effect-aware SEM returned no labelling")
    labels$SEM <- sem$labelling

    for (method in c("naive", "SEM")) {
      lr <- label_recovery_by_stratum(labels[[method]], treatment_time)
      label_rows[[length(label_rows) + 1L]] <-
        label_row_to_recovery(lr, method, rep_id)
    }

    fits <- list(
      naive = fit_effect_modified_hawkes(
        labels$naive, partition, design$grid, design$z_obs, post_window,
        treated_params, h_init = 0,
        maxit = fit_maxit
      ),
      SEM = sem$h_fit %||% failed_fit("Effect-aware SEM returned no h fit")
    )
    regimes <- list(
      z_flip_plus = design$z_flip_plus,
      z_flip_minus = design$z_flip_minus,
      z_all = design$z_all,
      z_none = design$z_none
    )
    truth_regimes <- simulate_regime_counts(
      partition, omega, post_window, regimes,
      control_params, treated_params, pre, design$grid, h_true,
      n_sim = n_truth,
      seed = structured_stage_seed(base_seed, 22L, rep_id, round(100 * h_true))
    )
    truth_targets <- contrast_effect_regimes(truth_regimes)
    names(truth_targets)[names(truth_targets) == "estimate"] <- "truth"
    for (method in names(fits)) {
      fit <- fits[[method]]
      if (!isTRUE(fit$converged)) {
        est <- data.frame(
          target = truth_targets$target,
          estimate = NA_real_,
          mc_se = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        est_regimes <- simulate_regime_counts(
          partition, omega, post_window, regimes,
          control_params,
          as.list(fit$par[c("mu", "alpha", "beta", "K")]),
          pre, design$grid, fit$par[["h"]],
          n_sim = n_cf,
          seed = structured_stage_seed(
            base_seed, 23L, rep_id,
            round(100 * h_true) + match(method, names(fits))
          )
        )
        est <- contrast_effect_regimes(est_regimes)
      }
      row <- left_join(est, truth_targets[, c("target", "truth")], by = "target")
      contrast_family <- c(
        psi_global = "global_1_0",
        psi_flip_plus = "single_cell_flip_X_plus",
        psi_flip_minus = "single_cell_flip_X_minus"
      )[row$target]
      causal_rows[[length(causal_rows) + 1L]] <- data.frame(
        method = method,
        contrast_family = unname(contrast_family),
        target = row$target,
        psi_estimate = row$estimate,
        psi_truth = row$truth,
        sim_id = rep_id,
        stringsAsFactors = FALSE
      )
    }
    list(
      ok = TRUE,
      labels = bind_rows(label_rows),
      psi = bind_rows(causal_rows),
      fit_status = bind_rows(lapply(names(fits), function(method) {
        fit_status_row(fits[[method]], method, rep_id)
      }))
    )
    }, error = function(e) {
      list(ok = FALSE, rep_id = rep_id, error = conditionMessage(e))
    })
  }
  parts <- collect_replication_parts(
    structured_lapply(seq_len(n_rep), run_effect_replication),
    sprintf("effect h=%.3f", h_true)
  )
  label_recovery_prefilter <- parts$labels
  psi_rows <- parts$psi
  fit_status <- parts$fit_status
  assert_complete_rows(
    label_recovery_prefilter, psi_rows, fit_status,
    c("global_1_0", "single_cell_flip_X_plus", "single_cell_flip_X_minus")
  )
} else {
  design <- make_geometry_transport_design(
    partition, path_seed = base_seed + 201L, observed_seed = base_seed + 202L
  )
  if (!is.finite(geometry_m) || !(geometry_m %in% design$m)) {
    stop(sprintf("PP_GEOMETRY_M=%s not in design path", as.character(geometry_m)))
  }
  alloc_idx <- match(geometry_m, design$m)
  z_est <- as.logical(design$allocations[[alloc_idx]])
  coarseness_val <- as.numeric(design$coarseness[[alloc_idx]])
  eval_regimes <- list(
    z_all = rep(TRUE, partition$n),
    z_none = rep(FALSE, partition$n)
  )
  run_geometry_replication <- function(rep_id) {
    tryCatch({
    label_rows <- list()
    causal_rows <- list()
    message(sprintf(
      "[structured-scenario] geometry m=%d C=%.3f rep %d/%d",
      geometry_m, coarseness_val, rep_id, n_rep
    ))
    pre <- simulate_pre(rep_id + 50000L + 1000L * alloc_idx)
    set.seed(structured_stage_seed(base_seed, 31L, rep_id, alloc_idx))
    post <- simulate_structured_catalogue(
      partition, omega, post_window, z_est,
      control_params, treated_params, pre, design$grid, h = 0
    )
    catalogue <- bind_catalogue_local(pre, post, design$grid, z_est, 0)
    labels <- make_oracle_and_naive_labels(catalogue, treatment_time)
    sem <- run_homogeneous_sem(
      catalogue, partition, z_est, omega, post_window,
      treatment_time, control_params, treated_params, sem_control,
      seed = structured_stage_seed(base_seed, 32L, rep_id, alloc_idx)
    )
    if (is.null(sem$labelling)) stop("Homogeneous SEM returned no labelling")
    labels$SEM <- sem$labelling

    for (method in c("naive", "SEM")) {
      lr <- label_recovery_by_stratum(labels[[method]], treatment_time)
      label_rows[[length(label_rows) + 1L]] <-
        label_row_to_recovery(lr, method, rep_id)
    }

    sem_par <- vapply(
      c("mu", "alpha", "beta", "K"),
      function(nm) as.numeric(sem$treated_params[[nm]] %||% NA_real_),
      numeric(1L)
    )
    fits <- list(
      naive = fit_effect_modified_hawkes(
        labels$naive, partition, design$grid, z_est,
        post_window, treated_params, homogeneous = TRUE, maxit = fit_maxit
      ),
      SEM = list(
        par = c(sem_par, h = 0),
        converged = all(is.finite(sem_par)) &&
          sem_par[["K"]] >= 0 && sem_par[["K"]] < 1,
        convergence = if (all(is.finite(sem_par))) 0L else 99L,
        message = ""
      )
    )

    cells <- design$focal_cells
    truth_regimes <- simulate_regime_counts(
      partition, omega, post_window, eval_regimes,
      control_params, treated_params, pre, design$grid, h = 0,
      region_cells = cells, n_sim = n_truth,
      seed = structured_stage_seed(base_seed, 33L, rep_id, alloc_idx)
    )
    zero <- truth_regimes$mean[truth_regimes$allocation == "z_none"]
    all_t <- truth_regimes$mean[truth_regimes$allocation == "z_all"]
    psi_truth <- as.numeric(all_t - zero)
    for (method in names(fits)) {
      fit <- fits[[method]]
      if (!isTRUE(fit$converged)) {
        psi_hat <- NA_real_
      } else {
        regimes <- simulate_regime_counts(
          partition, omega, post_window, eval_regimes,
          control_params, as.list(fit$par[c("mu", "alpha", "beta", "K")]),
          pre, design$grid, h = 0, region_cells = cells, n_sim = n_cf,
          seed = structured_stage_seed(
            base_seed, 34L, rep_id,
            100L * alloc_idx + match(method, names(fits))
          )
        )
        psi_hat <- as.numeric(
          regimes$mean[regimes$allocation == "z_all"] -
            regimes$mean[regimes$allocation == "z_none"]
        )
      }
      causal_rows[[length(causal_rows) + 1L]] <- data.frame(
        method = method,
        contrast_family = "focal_all_1_0",
        target = "psi_focal_all_1_0",
        psi_estimate = psi_hat,
        psi_truth = psi_truth,
        sim_id = rep_id,
        stringsAsFactors = FALSE
      )
    }
    list(
      ok = TRUE,
      labels = bind_rows(label_rows),
      psi = bind_rows(causal_rows),
      fit_status = bind_rows(lapply(names(fits), function(method) {
        fit_status_row(fits[[method]], method, rep_id)
      }))
    )
    }, error = function(e) {
      list(ok = FALSE, rep_id = rep_id, error = conditionMessage(e))
    })
  }
  parts <- collect_replication_parts(
    structured_lapply(seq_len(n_rep), run_geometry_replication),
    sprintf("geometry m=%d", geometry_m)
  )
  label_recovery_prefilter <- parts$labels
  psi_rows <- parts$psi
  fit_status <- parts$fit_status
  assert_complete_rows(
    label_recovery_prefilter, psi_rows, fit_status,
    "focal_all_1_0"
  )
}

lr_df <- label_recovery_prefilter
psi_df <- psi_rows
class_metrics <- summarize_class(lr_df)
support_contrast_summary <- summarize_support(psi_df)

out_root <- PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir)
rds_path <- file.path(out_root, paste0(run_basename, ".rds"))
dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)

result <- list(
  class_metrics = class_metrics,
  label_recovery_prefilter = lr_df,
  fit_status = fit_status,
  support_contrast_summary = support_contrast_summary,
  support_contrast_df = psi_df,
  all_nothing_ATE = if (!is.null(psi_df) && nrow(psi_df) > 0L) {
    mean(psi_df$psi_truth, na.rm = TRUE)
  } else NA_real_,
  decay_validation = NULL,
  summary_df = NULL,
  config = list(
    SCENARIO_ID = scenario_id,
    SCENARIO_FAMILY = family,
    N_SIMS = n_rep,
    REP_WORKERS = rep_workers,
    TRUTH_SIMS = n_truth,
    FORWARD_SIMS = n_cf,
    TRUE_MU = true_mu,
    TARGET_POINTS = target_points,
    H_TRUE = h_true,
    GEOMETRY_M = geometry_m,
    COARSENESS = coarseness_val,
    hawkes_par_1 = control_params,
    hawkes_par_2 = treated_params,
    NX = 10L, NY = 10L, TREAT_PROP = 0.5,
    structured_scenario = TRUE
  )
)
saveRDS(result, rds_path)
message("[structured-scenario] wrote ", rds_path)
if (!is.null(class_metrics)) print(class_metrics)
if (!is.null(support_contrast_summary)) print(support_contrast_summary)
