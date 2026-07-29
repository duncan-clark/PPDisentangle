#!/usr/bin/env Rscript
# Structured robustness studies:
#   1. Binary-covariate treatment-effect modification.
#   2. Transport across balanced allocation geometries.
#
# Examples:
#   Rscript inst/sim_study/sim_study_structured_robustness.R --test
#   Rscript inst/sim_study/sim_study_structured_robustness.R --study effect --pilot-only
#   Rscript inst/sim_study/sim_study_structured_robustness.R --sims 100

suppressPackageStartupMessages({
  library(PPDisentangle)
  library(spatstat)
  library(dplyr)
  library(ggplot2)
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
arg_num <- function(flag, default) {
  x <- suppressWarnings(as.numeric(get_arg(flag, default)))
  if (!is.finite(x) || is.na(x)) default else x
}

script_dir <- {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1L]])))
  else normalizePath(file.path(getwd(), "inst", "sim_study"), mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(script_dir, "structured_robustness_utils.R"), local = FALSE)

test_mode <- "--test" %in% args
pilot_only <- "--pilot-only" %in% args
skip_pilot_gate <- "--skip-pilot-gate" %in% args
study <- tolower(get_arg("--study", "both"))
if (!study %in% c("both", "effect", "geometry")) {
  stop("--study must be one of: both, effect, geometry")
}

base_seed <- arg_int("--seed", 20260714L)
n_rep <- arg_int("--sims", if (test_mode) 1L else 32L)
pilot_reps <- arg_int("--pilot-reps", if (test_mode) 1L else 20L)
n_truth <- arg_int("--truth-sims", if (test_mode) 2L else 500L)
n_cf <- arg_int("--forward-sims", if (test_mode) 2L else 100L)
sem_inner <- arg_int("--sem-inner", if (test_mode) 5L else 1000L)
sem_outer <- arg_int("--sem-outer", if (test_mode) 1L else 3L)
sem_labels <- arg_int("--sem-labellings", if (test_mode) 1L else 10L, 0L)
sem_props <- arg_int("--sem-props", if (test_mode) 1L else 20L)
# One effect-aware SEM pass per replication (no homogeneous comparator SEM).
h_rounds <- arg_int("--h-sem-rounds", 1L)
fit_maxit <- arg_int("--fit-maxit", if (test_mode) 30L else 1000L)
n_workers <- arg_int(
  "--workers",
  {
    slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
    env <- suppressWarnings(as.integer(Sys.getenv("PP_STRUCTURED_WORKERS", "")))
    if (is.finite(slurm) && !is.na(slurm) && slurm >= 1L) slurm
    else if (is.finite(env) && !is.na(env) && env >= 1L) env
    else max(1L, parallel::detectCores(logical = FALSE) - 1L)
  }
)
if (test_mode) n_workers <- 1L

structured_lapply <- function(X, FUN) {
  n <- length(X)
  if (n < 1L) return(list())
  cores <- max(1L, min(as.integer(n_workers), n))
  if (cores <= 1L) return(lapply(X, FUN))
  message(sprintf("[structured] parallel workers=%d over %d jobs", cores, n))
  parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
}

omega <- spatstat.geom::owin(c(0, 100), c(0, 100))
partition <- spatstat.geom::quadrats(X = omega, nx = 10, ny = 10)
treatment_time <- 10
end_time <- 110
post_window <- c(treatment_time, end_time)
control_K <- 0.8
treated_K <- 0.2
target_points <- arg_num("--target-points", 2500)
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
  param_update_cadence = 10L,
  proposal_update_cadence = 1L,
  update_control_params = FALSE,
  include_starting_data = TRUE,
  update_starting_data = TRUE,
  maxit = fit_maxit,
  verbose = FALSE
)

out_root <- PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir)
robustness_generated_dir <- file.path(out_root, "generated", "robustness")
generated_dir <- file.path(robustness_generated_dir, "structured")
figure_dir <- file.path(robustness_generated_dir, "figures")
table_dir <- file.path(robustness_generated_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
run_id <- get_arg("--output", paste0("structured_", format(Sys.time(), "%Y%m%d_%H%M%S")))

message(sprintf(
  "[structured] study=%s reps=%d pilot=%d truth_sims=%d forward_sims=%d workers=%d h_sem_rounds=%d",
  study, n_rep, pilot_reps, n_truth, n_cf, n_workers, h_rounds
))
message(sprintf(
  "[structured] conditioning: truths and estimates share each replication's pre-treatment history"
))
message(paste0(
  "[structured] design: effect = 1 SEM + naive per (h,rep); ",
  "geometry = simulate/SEM/fit under each coarseness regime"
))

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

bind_catalogue <- function(pre, post, partition, grid, z_obs, h) {
  pre <- attach_design_columns(pre, partition, grid, z_obs, h = 0)
  pre$process <- "control"
  post <- attach_design_columns(post, partition, grid, z_obs, h = h)
  common <- union(names(pre), names(post))
  for (nm in setdiff(common, names(pre))) pre[[nm]] <- NA
  for (nm in setdiff(common, names(post))) post[[nm]] <- NA
  out <- rbind(pre[, common, drop = FALSE], post[, common, drop = FALSE])
  out <- out[order(out$t), , drop = FALSE]
  # adaptive_SEM treats an existing inferred_process column as its starting
  # labelling, so it must be complete rather than NA for post-treatment rows.
  out$inferred_process <- ifelse(
    out$t < treatment_time, "control", as.character(out$location_process)
  )
  out
}

fit_rows_from_list <- function(fits, rep_id, h_true) {
  do.call(rbind, lapply(names(fits), function(key) {
    fit <- fits[[key]]
    bits <- strsplit(key, "__", fixed = TRUE)[[1L]]
    hhat <- fit$par[["h"]]
    hse <- fit$h_se %||% NA_real_
    data.frame(
      replication = rep_id,
      h_true = h_true,
      model = bits[[1L]],
      method = bits[[2L]],
      h_hat = hhat,
      h_se = hse,
      covered = is.finite(hhat) && is.finite(hse) &&
        h_true >= hhat - 1.96 * hse && h_true <= hhat + 1.96 * hse,
      mu1_hat = fit$par[["mu"]],
      alpha1_hat = fit$par[["alpha"]],
      beta1_hat = fit$par[["beta"]],
      K1_hat = fit$par[["K"]],
      converged = isTRUE(fit$converged),
      hessian_ok = isTRUE(fit$hessian_ok),
      boundary_hit = is.finite(hhat) && abs(hhat) > 10,
      stringsAsFactors = FALSE
    )
  }))
}

summarize_h <- function(h_df) {
  h_df |>
    group_by(.data$h_true, .data$model, .data$method) |>
    summarise(
      n_rep = n(),
      n_ok = sum(.data$converged & is.finite(.data$h_hat)),
      mean_estimate = mean(.data$h_hat[.data$converged], na.rm = TRUE),
      bias = mean(.data$h_hat[.data$converged] - .data$h_true[.data$converged], na.rm = TRUE),
      rmse = sqrt(mean((.data$h_hat[.data$converged] -
        .data$h_true[.data$converged])^2, na.rm = TRUE)),
      mean_hessian_se = mean(.data$h_se[.data$converged], na.rm = TRUE),
      wald_coverage = mean(.data$covered[.data$converged], na.rm = TRUE),
      convergence_rate = mean(.data$converged),
      boundary_hit_rate = mean(.data$boundary_hit, na.rm = TRUE),
      corr_h_mu1 = if (sum(.data$converged) > 2L) {
        suppressWarnings(cor(.data$h_hat[.data$converged],
                             .data$mu1_hat[.data$converged]))
      } else NA_real_,
      .groups = "drop"
    )
}

contrast_effect_regimes <- function(regime_df) {
  val <- setNames(regime_df$mean, regime_df$allocation)
  se <- setNames(regime_df$mc_se, regime_df$allocation)
  need <- c("z_flip_plus", "z_flip_minus", "z_all", "z_none")
  if (!all(need %in% names(val))) {
    stop(
      "contrast_effect_regimes expects allocations: ",
      paste(need, collapse = ", ")
    )
  }
  targets <- c(
    psi_global = val[["z_all"]] - val[["z_none"]],
    psi_flip_plus = val[["z_flip_plus"]] - val[["z_none"]],
    psi_flip_minus = val[["z_flip_minus"]] - val[["z_none"]]
  )
  target_se <- c(
    psi_global = sqrt(se[["z_all"]]^2 + se[["z_none"]]^2),
    psi_flip_plus = sqrt(se[["z_flip_plus"]]^2 + se[["z_none"]]^2),
    psi_flip_minus = sqrt(se[["z_flip_minus"]]^2 + se[["z_none"]]^2)
  )
  data.frame(target = names(targets), estimate = unname(targets),
             mc_se = unname(target_se), stringsAsFactors = FALSE)
}

run_oracle_pilot <- function(design) {
  h_values <- c(0, 0.3, 0.6)
  h_rows <- list()
  profile_rows <- list()
  count_rows <- list()
  delta_rows <- list()
  k <- 0L
  for (h_true in h_values) {
    for (rep_id in seq_len(pilot_reps)) {
      pre <- simulate_pre(rep_id + 10000L * match(h_true, h_values))
      set.seed(structured_stage_seed(base_seed, 11L, rep_id, round(100 * h_true)))
      post <- simulate_structured_catalogue(
        partition, omega, post_window, design$z_obs,
        control_params, treated_params, pre, design$grid, h_true
      )
      catalogue <- bind_catalogue(pre, post, partition, design$grid, design$z_obs, h_true)
      oracle <- make_oracle_and_naive_labels(catalogue, treatment_time)$oracle
      fit <- fit_effect_modified_hawkes(
        oracle, partition, design$grid, design$z_obs, post_window,
        treated_params, h_init = 0, maxit = fit_maxit
      )
      k <- k + 1L
      h_rows[[k]] <- fit_rows_from_list(
        list(heterogeneous__oracle = fit), rep_id, h_true
      )
      post_treated <- oracle[
        oracle$t >= treatment_time & oracle$process == "treated", , drop = FALSE
      ]
      count_rows[[k]] <- data.frame(
        replication = rep_id, h_true = h_true,
        n_X_plus = sum(post_treated$X == 1),
        n_X_minus = sum(post_treated$X == -1)
      )
      regimes <- list(
        z_plus = design$z_plus, z_minus = design$z_minus,
        z_none = design$z_none
      )
      truth_regimes <- simulate_regime_counts(
        partition, omega, post_window, regimes,
        control_params, treated_params, pre, design$grid, h_true,
        n_sim = n_truth,
        seed = structured_stage_seed(base_seed, 12L, rep_id, round(100 * h_true))
      )
      if (isTRUE(fit$converged)) {
        estimate_regimes <- simulate_regime_counts(
          partition, omega, post_window, regimes,
          control_params, as.list(fit$par[c("mu", "alpha", "beta", "K")]),
          pre, design$grid, fit$par[["h"]],
          n_sim = n_cf,
          seed = structured_stage_seed(base_seed, 13L, rep_id, round(100 * h_true))
        )
        delta_hat <- with(
          estimate_regimes,
          mean[allocation == "z_plus"] - mean[allocation == "z_minus"]
        )
      } else {
        delta_hat <- NA_real_
      }
      delta_truth <- with(
        truth_regimes,
        mean[allocation == "z_plus"] - mean[allocation == "z_minus"]
      )
      delta_rows[[k]] <- data.frame(
        replication = rep_id, h_true = h_true,
        estimate = delta_hat, truth = delta_truth,
        mc_se = sqrt(sum(truth_regimes$mc_se[
          truth_regimes$allocation %in% c("z_plus", "z_minus")
        ]^2))
      )
      if (rep_id == 1L) {
        profile <- profile_effect_h(
          fit, oracle, partition, design$grid, design$z_obs, post_window,
          h_grid = seq(-0.5, 1.1, length.out = if (test_mode) 3L else 25L),
          maxit = fit_maxit
        )
        profile$h_true <- h_true
        profile_rows[[length(profile_rows) + 1L]] <- profile
      }
    }
  }
  h_df <- bind_rows(h_rows)
  summary <- summarize_h(h_df)
  list(
    h = h_df,
    summary = summary,
    profiles = bind_rows(profile_rows),
    event_counts = bind_rows(count_rows),
    delta_het = bind_rows(delta_rows),
    delta_het_summary = causal_error_summary(
      bind_rows(delta_rows), group_cols = "h_true"
    )
  )
}

run_effect_study <- function() {
  design <- make_effect_modification_design(partition, seed = base_seed + 101L)
  pilot <- run_oracle_pilot(design)
  pilot_csv <- file.path(table_dir, "effect_modification_oracle_pilot.csv")
  write.csv(pilot$summary, pilot_csv, row.names = FALSE)
  gate <- pilot$summary |>
    filter(.data$model == "heterogeneous", .data$method == "oracle") |>
    summarise(
      pass = all(.data$convergence_rate >= 0.8) &&
        all(abs(.data$bias) <= 0.25) && all(.data$rmse <= 0.4)
    ) |>
    pull(.data$pass)
  if (test_mode) gate <- TRUE
  if (!isTRUE(gate) && !skip_pilot_gate) {
    stop(
      "Oracle pilot gate failed (requires convergence >= .8, |bias| <= .25, RMSE <= .4). ",
      "Inspect ", pilot_csv, " or use --skip-pilot-gate after review."
    )
  }
  if (pilot_only) {
    return(list(design = design, pilot = pilot, pilot_only = TRUE))
  }

  h_values <- c(0, 0.3, 0.6)
  jobs <- expand.grid(
    h_true = h_values,
    rep_id = seq_len(n_rep),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "[effect] %d data realizations (%d h x %d reps); 1 SEM + naive fit each; workers=%d",
    nrow(jobs), length(h_values), n_rep, n_workers
  ))

  run_one_effect <- function(job_i) {
    h_true <- jobs$h_true[[job_i]]
    rep_id <- as.integer(jobs$rep_id[[job_i]])
    message(sprintf("[effect] h=%.1f replication %d/%d", h_true, rep_id, n_rep))
    pre <- simulate_pre(rep_id + 1000L * match(h_true, h_values))
    set.seed(structured_stage_seed(base_seed, 20L, rep_id, round(100 * h_true)))
    post <- simulate_structured_catalogue(
      partition, omega, post_window, design$z_obs,
      control_params, treated_params, pre, design$grid, h_true
    )
    catalogue <- bind_catalogue(pre, post, partition, design$grid, design$z_obs, h_true)
    labels <- make_oracle_and_naive_labels(catalogue, treatment_time)
    # Single SEM labelling per realization (one effect-aware pass).
    sem <- tryCatch(
      run_effect_aware_sem(
        catalogue, partition, design$grid, design$z_obs, omega,
        post_window, treatment_time, control_params, treated_params,
        h_init = 0, interaction_rounds = h_rounds,
        sem_control = sem_control,
        seed = structured_stage_seed(base_seed, 21L, rep_id, round(100 * h_true))
      ),
      error = function(e) {
        warning(sprintf(
          "Effect-aware SEM failed for h=%.1f replication %d: %s",
          h_true, rep_id, conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(sem)) labels$SEM <- sem$labelling

    fits <- list()
    for (method in intersect(c("oracle", "naive", "SEM"), names(labels))) {
      fits[[paste0("heterogeneous__", method)]] <- fit_effect_modified_hawkes(
        labels[[method]], partition, design$grid, design$z_obs, post_window,
        treated_params, h_init = if (method == "SEM" && !is.null(sem)) sem$h else 0,
        maxit = fit_maxit
      )
    }
    h_row <- fit_rows_from_list(fits, rep_id, h_true)

    label_rows <- list()
    for (method in names(labels)) {
      lr <- label_recovery_by_stratum(labels[[method]], treatment_time)
      lr$replication <- rep_id
      lr$h_true <- h_true
      lr$method <- method
      label_rows[[length(label_rows) + 1L]] <- lr
    }

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
    causal_rows <- list()
    for (fit_name in names(fits)) {
      fit <- fits[[fit_name]]
      bits <- strsplit(fit_name, "__", fixed = TRUE)[[1L]]
      if (!isTRUE(fit$converged)) {
        est_targets <- data.frame(
          target = truth_targets$target, estimate = NA_real_, mc_se = NA_real_
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
            round(100 * h_true) + match(fit_name, names(fits))
          )
        )
        est_targets <- contrast_effect_regimes(est_regimes)
      }
      row <- left_join(est_targets, truth_targets[, c("target", "truth")], by = "target")
      row$replication <- rep_id
      row$h_true <- h_true
      row$model <- bits[[1L]]
      row$method <- bits[[2L]]
      causal_rows[[length(causal_rows) + 1L]] <- row
    }
    list(
      h = h_row,
      causal = bind_rows(causal_rows),
      labels = bind_rows(label_rows)
    )
  }

  effect_parts <- structured_lapply(seq_len(nrow(jobs)), run_one_effect)
  h_df <- bind_rows(lapply(effect_parts, `[[`, "h"))
  causal_df <- bind_rows(lapply(effect_parts, `[[`, "causal"))
  label_df <- bind_rows(lapply(effect_parts, `[[`, "labels"))
  h_summary <- summarize_h(h_df)
  causal_summary <- causal_error_summary(
    causal_df, group_cols = c("h_true", "model", "method", "target")
  )
  label_summary <- label_df |>
    group_by(.data$h_true, .data$method, .data$stratum) |>
    summarise(
      accuracy = mean(.data$accuracy, na.rm = TRUE),
      balanced_accuracy = mean(.data$balanced_accuracy, na.rm = TRUE),
      recall_control = mean(.data$recall_control, na.rm = TRUE),
      recall_treated = mean(.data$recall_treated, na.rm = TRUE),
      .groups = "drop"
    )
  list(
    design = design, pilot = pilot, h = h_df, h_summary = h_summary,
    causal = causal_df, causal_summary = causal_summary,
    labels = label_df, label_summary = label_summary,
    pilot_only = FALSE
  )
}

run_geometry_truth_pilot <- function(design) {
  pre <- simulate_pre(90001L)
  regions <- list(focal = design$focal_cells, domain = seq_len(partition$n))
  out <- list()
  for (region in names(regions)) {
    truth <- simulate_regime_counts(
      partition, omega, post_window,
      c(design$allocations, list(z_none = rep(FALSE, partition$n))),
      control_params, treated_params, pre, design$grid, h = 0,
      region_cells = regions[[region]], n_sim = n_truth,
      seed = structured_stage_seed(base_seed, 30L, match(region, names(regions)))
    )
    zero <- truth$mean[truth$allocation == "z_none"]
    zero_se <- truth$mc_se[truth$allocation == "z_none"]
    truth <- truth[truth$allocation != "z_none", , drop = FALSE]
    truth$contrast <- truth$mean - zero
    truth$contrast_mc_se <- sqrt(truth$mc_se^2 + zero_se^2)
    truth$region <- region
    truth$m <- design$m
    truth$coarseness <- unname(design$coarseness)
    out[[region]] <- truth
  }
  pilot <- bind_rows(out)
  focal <- pilot[pilot$region == "focal", ]
  domain <- pilot[pilot$region == "domain", ]
  focal_range <- diff(range(focal$contrast))
  domain_range <- diff(range(domain$contrast))
  isolated_share <- max(abs(diff(focal$contrast))) / max(focal_range, 1e-12)
  diagnostics <- data.frame(
    focal_range = focal_range,
    max_focal_mc_se = max(focal$contrast_mc_se),
    domain_range = domain_range,
    max_domain_mc_se = max(domain$contrast_mc_se),
    largest_step_share = isolated_share,
    coarseness_min = min(design$coarseness),
    coarseness_max = max(design$coarseness),
    pass = focal_range > 2 * max(focal$contrast_mc_se) &&
      min(design$coarseness) == 0 && max(design$coarseness) == 1
  )
  list(truth = pilot, diagnostics = diagnostics, pre_history = pre)
}

geometry_global_contrast <- function(regime_df, m, coarseness, region) {
  zero <- regime_df$mean[regime_df$allocation == "z_none"]
  zero_se <- regime_df$mc_se[regime_df$allocation == "z_none"]
  all_treated <- regime_df$mean[regime_df$allocation == "z_all"]
  all_se <- regime_df$mc_se[regime_df$allocation == "z_all"]
  data.frame(
    target = "psi_global",
    m = as.integer(m),
    coarseness = as.numeric(coarseness),
    region = region,
    estimate = as.numeric(all_treated - zero),
    mc_se = sqrt(as.numeric(all_se)^2 + as.numeric(zero_se)^2),
    stringsAsFactors = FALSE
  )
}

run_geometry_study <- function() {
  design <- make_geometry_transport_design(
    partition, path_seed = base_seed + 201L, observed_seed = base_seed + 202L
  )
  pilot <- run_geometry_truth_pilot(design)
  write.csv(pilot$truth, file.path(table_dir, "geometry_truth_pilot.csv"), row.names = FALSE)
  write.csv(pilot$diagnostics, file.path(table_dir, "geometry_truth_pilot_diagnostics.csv"),
            row.names = FALSE)
  if (!isTRUE(pilot$diagnostics$pass) && !skip_pilot_gate && !test_mode) {
    stop(
      "Geometry truth pilot gate failed: focal variation is not separated from MC error. ",
      "Inspect geometry_truth_pilot_diagnostics.csv or use --skip-pilot-gate after review."
    )
  }
  if (pilot_only) {
    return(list(design = design, pilot = pilot, pilot_only = TRUE))
  }

  n_alloc <- length(design$allocations)
  jobs <- expand.grid(
    alloc_idx = seq_len(n_alloc),
    rep_id = seq_len(n_rep),
    KEEP.OUT.ATTRS = FALSE
  )
  message(sprintf(
    paste0(
      "[geometry] %d jobs (%d coarseness regimes x %d reps); ",
      "simulate/SEM/fit under each regime; workers=%d"
    ),
    nrow(jobs), n_alloc, n_rep, n_workers
  ))

  # Evaluation target is fixed (all-treated vs all-control). Estimation design
  # varies: each job uses one path allocation as the observed treatment regime.
  eval_regimes <- list(
    z_all = rep(TRUE, partition$n),
    z_none = rep(FALSE, partition$n)
  )

  run_one_geometry <- function(job_i) {
    alloc_idx <- as.integer(jobs$alloc_idx[[job_i]])
    rep_id <- as.integer(jobs$rep_id[[job_i]])
    z_est <- as.logical(design$allocations[[alloc_idx]])
    m_val <- as.integer(design$m[[alloc_idx]])
    c_val <- as.numeric(design$coarseness[[alloc_idx]])
    message(sprintf(
      "[geometry] m=%d C=%.2f replication %d/%d",
      m_val, c_val, rep_id, n_rep
    ))

    pre <- simulate_pre(rep_id + 50000L + 1000L * alloc_idx)
    set.seed(structured_stage_seed(base_seed, 31L, rep_id, alloc_idx))
    post <- simulate_structured_catalogue(
      partition, omega, post_window, z_est,
      control_params, treated_params, pre, design$grid, h = 0
    )
    catalogue <- bind_catalogue(pre, post, partition, design$grid, z_est, 0)
    labels <- make_oracle_and_naive_labels(catalogue, treatment_time)
    sem <- tryCatch(
      run_homogeneous_sem(
        catalogue, partition, z_est, omega, post_window,
        treatment_time, control_params, treated_params, sem_control,
        seed = structured_stage_seed(base_seed, 32L, rep_id, alloc_idx)
      ),
      error = function(e) {
        warning(sprintf(
          "Geometry SEM failed for m=%d replication %d: %s",
          m_val, rep_id, conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(sem)) labels$SEM <- sem$labelling
    fits <- lapply(names(labels), function(method) {
      fit_effect_modified_hawkes(
        labels[[method]], partition, design$grid, z_est,
        post_window, treated_params, homogeneous = TRUE, maxit = fit_maxit
      )
    })
    names(fits) <- names(labels)

    label_rows <- list()
    parameter_rows <- list()
    for (method in names(labels)) {
      lr <- label_recovery_by_stratum(labels[[method]], treatment_time)
      lr$replication <- rep_id
      lr$method <- method
      lr$m <- m_val
      lr$coarseness <- c_val
      label_rows[[length(label_rows) + 1L]] <- lr
      fit <- fits[[method]]
      parameter_rows[[length(parameter_rows) + 1L]] <- data.frame(
        replication = rep_id, method = method, m = m_val, coarseness = c_val,
        converged = isTRUE(fit$converged),
        mu1_hat = fit$par[["mu"]], alpha1_hat = fit$par[["alpha"]],
        beta1_hat = fit$par[["beta"]], K1_hat = fit$par[["K"]]
      )
    }

    causal_rows <- list()
    for (region in c("focal", "domain")) {
      cells <- if (region == "focal") design$focal_cells else seq_len(partition$n)
      truth_regimes <- simulate_regime_counts(
        partition, omega, post_window, eval_regimes,
        control_params, treated_params, pre, design$grid, h = 0,
        region_cells = cells, n_sim = n_truth,
        seed = structured_stage_seed(
          base_seed, 33L, rep_id,
          10L * alloc_idx + match(region, c("focal", "domain"))
        )
      )
      truth <- geometry_global_contrast(truth_regimes, m_val, c_val, region)
      for (method in names(fits)) {
        fit <- fits[[method]]
        if (!isTRUE(fit$converged)) {
          est <- truth
          est$estimate <- NA_real_
          est$mc_se <- NA_real_
        } else {
          regimes <- simulate_regime_counts(
            partition, omega, post_window, eval_regimes,
            control_params, as.list(fit$par[c("mu", "alpha", "beta", "K")]),
            pre, design$grid, h = 0, region_cells = cells, n_sim = n_cf,
            seed = structured_stage_seed(
              base_seed, 34L, rep_id,
              100L * alloc_idx + 10L * match(method, names(fits)) +
                match(region, c("focal", "domain"))
            )
          )
          est <- geometry_global_contrast(regimes, m_val, c_val, region)
        }
        est$truth <- truth$estimate
        est$replication <- rep_id
        est$method <- method
        causal_rows[[length(causal_rows) + 1L]] <- est
      }
    }
    list(
      causal = bind_rows(causal_rows),
      labels = bind_rows(label_rows),
      parameters = bind_rows(parameter_rows)
    )
  }

  geometry_parts <- structured_lapply(seq_len(nrow(jobs)), run_one_geometry)
  causal <- bind_rows(lapply(geometry_parts, `[[`, "causal"))
  labels <- bind_rows(lapply(geometry_parts, `[[`, "labels"))
  parameters <- bind_rows(lapply(geometry_parts, `[[`, "parameters"))
  summary <- causal_error_summary(
    causal, group_cols = c("region", "method", "m", "coarseness", "target")
  )
  label_summary <- labels |>
    filter(.data$stratum == "all") |>
    group_by(.data$m, .data$coarseness, .data$method) |>
    summarise(
      accuracy = mean(.data$accuracy, na.rm = TRUE),
      balanced_accuracy = mean(.data$balanced_accuracy, na.rm = TRUE),
      recall_control = mean(.data$recall_control, na.rm = TRUE),
      recall_treated = mean(.data$recall_treated, na.rm = TRUE),
      .groups = "drop"
    )
  list(
    design = design, pilot = pilot, causal = causal, causal_summary = summary,
    labels = labels, label_summary = label_summary,
    parameters = parameters, pilot_only = FALSE
  )
}

allocation_map_df <- function(grid, allocations, labels, focal_cells = integer()) {
  out <- bind_rows(lapply(seq_along(allocations), function(i) {
    transform(
      grid,
      allocation = labels[[i]],
      treated = ifelse(allocations[[i]], "Treated", "Control"),
      focal = cell %in% focal_cells
    )
  }))
  out$allocation <- factor(out$allocation, levels = labels)
  out
}

plot_allocation_maps <- function(map_df, nrow = 1L, show_X = FALSE) {
  fill_var <- if (show_X) "treated" else "treated"
  # Outline only the marked cell(s): geom_tile borders stroke across neighbors;
  # geom_rect with unit cell bounds keeps the stroke on that one tile.
  focal_df <- map_df[as.logical(map_df$focal), c("col", "row", "allocation"), drop = FALSE]
  focal_df <- unique(focal_df)
  # Keep the stroke inside the tile so neighbors are not painted black.
  inset <- 0.12
  ggplot(map_df, aes(x = .data$col, y = .data$row, fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 0.15, width = 1, height = 1) +
    {
      if (nrow(focal_df)) {
        geom_rect(
          data = focal_df,
          aes(
            xmin = .data$col - 0.5 + inset, xmax = .data$col + 0.5 - inset,
            ymin = .data$row - 0.5 + inset, ymax = .data$row + 0.5 - inset
          ),
          fill = NA, color = "black", linewidth = 0.35,
          inherit.aes = FALSE
        )
      } else {
        NULL
      }
    } +
    facet_wrap(vars(.data$allocation), nrow = nrow) +
    scale_fill_manual(values = c(
      Control = "#D9E2EC", Treated = "#D1495B",
      `X = -1` = "#74A9CF", `X = +1` = "#FCAE91"
    )) +
    coord_equal(expand = FALSE) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(), strip.text = element_text(size = 13),
          legend.position = "bottom")
}

save_ggplot_grid <- function(plots, file, widths, heights, layout_matrix,
                             width, height, dpi = 300) {
  if (grepl("[.]pdf$", file, ignore.case = TRUE)) {
    grDevices::pdf(file, width = width, height = height, bg = "white",
                   useDingbats = FALSE)
  } else {
    grDevices::png(file, width = width, height = height, units = "in",
                   res = dpi, bg = "white")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  lay <- grid::grid.layout(
    nrow(layout_matrix), ncol(layout_matrix),
    widths = grid::unit(widths, "null"), heights = grid::unit(heights, "null")
  )
  vp <- grid::viewport(layout = lay)
  grid::pushViewport(vp)
  for (i in seq_along(plots)) {
    pos <- which(layout_matrix == i, arr.ind = TRUE)
    print(plots[[i]], vp = grid::viewport(
      layout.pos.row = range(pos[, "row"]),
      layout.pos.col = range(pos[, "col"])
    ), newpage = FALSE)
  }
  grid::popViewport()
  invisible(file)
}

make_effect_figure <- function(res) {
  d <- res$design
  x_map <- transform(
    d$grid, allocation = "Covariate X",
    treated = ifelse(X == 1, "X = +1", "X = -1"),
    focal = cell %in% c(d$flip_plus_cell, d$flip_minus_cell)
  )
  z_map <- allocation_map_df(
    d$grid, list(d$z_obs), c("Observed z[obs]"),
    focal_cells = c(d$flip_plus_cell, d$flip_minus_cell)
  )
  maps <- bind_rows(x_map, z_map)
  maps$allocation <- factor(
    as.character(maps$allocation),
    levels = c("Covariate X", "Observed z[obs]")
  )
  pA <- plot_allocation_maps(maps, nrow = 1L) +
    labs(title = sprintf(
      "A  Fixed X and z[obs] (outlined single-cell flips: X=+1 cell %d, X=-1 cell %d)",
      d$flip_plus_cell, d$flip_minus_cell
    ))

  bias_df <- res$causal_summary |>
    filter(
      .data$model == "heterogeneous",
      .data$method %in% c("naive", "SEM"),
      .data$target %in% c("psi_global", "psi_flip_plus", "psi_flip_minus")
    ) |>
    mutate(
      method = factor(.data$method, levels = c("naive", "SEM")),
      target = factor(
        .data$target,
        levels = c("psi_global", "psi_flip_plus", "psi_flip_minus"),
        labels = c(
          "All-or-nothing DAITE",
          "Single X=+1 cell flip",
          "Single X=-1 cell flip"
        )
      ),
      se = ifelse(is.finite(.data$mcse_bias), .data$mcse_bias, 0)
    )
  pB <- ggplot(bias_df, aes(.data$h_true, .data$bias,
                            color = .data$method, shape = .data$method)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_linerange(aes(ymin = .data$bias - .data$se,
                       ymax = .data$bias + .data$se),
                   position = position_dodge(width = 0.04)) +
    geom_line(position = position_dodge(width = 0.04)) +
    geom_point(position = position_dodge(width = 0.04), size = 2) +
    facet_wrap(vars(.data$target), scales = "free_y", nrow = 1L) +
    scale_color_manual(values = c(naive = "#E07A5F", SEM = "#3D405B")) +
    scale_shape_manual(values = c(naive = 16, SEM = 17)) +
    labs(
      title = "B  Estimated DAITE bias",
      x = "True h",
      y = "DAITE bias (estimate - truth)",
      color = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 13) + theme(legend.position = "bottom")

  list(pA, pB)
}

make_geometry_figure <- function(res) {
  d <- res$design
  labels <- sprintf("m=%d; C=%.2f", d$m, d$coarseness)
  maps <- allocation_map_df(d$grid, d$allocations, labels, d$focal_cells)
  pA <- plot_allocation_maps(maps, nrow = 1L) +
    labs(title = "A  Estimation regimes along the coarseness path (each used as z for simulate/SEM/fit)")
  bias_df <- res$causal_summary |>
    filter(
      .data$region == "focal",
      .data$method %in% c("naive", "SEM"),
      is.na(.data$target) | .data$target == "psi_global"
    ) |>
    mutate(
      method = factor(.data$method, levels = c("naive", "SEM")),
      se = ifelse(is.finite(.data$mcse_bias), .data$mcse_bias, 0)
    )
  pB <- ggplot(bias_df, aes(.data$coarseness, .data$bias,
                            color = .data$method, shape = .data$method)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_linerange(aes(ymin = .data$bias - .data$se,
                       ymax = .data$bias + .data$se),
                   position = position_dodge(width = 0.02)) +
    geom_line(position = position_dodge(width = 0.02)) +
    geom_point(position = position_dodge(width = 0.02), size = 2) +
    scale_color_manual(values = c(naive = "#E07A5F", SEM = "#3D405B")) +
    scale_shape_manual(values = c(naive = 16, SEM = 17)) +
    labs(
      title = "B  Estimated all-or-nothing DAITE bias by estimation-regime coarseness",
      x = "Coarseness C(z) of the estimation regime",
      y = "DAITE bias (estimate - truth)",
      color = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 13) + theme(legend.position = "bottom")
  list(pA, pB)
}

write_compact_tex <- function(df, path, caption, label, digits = 3L) {
  num <- vapply(df, is.numeric, logical(1L))
  df[num] <- lapply(df[num], function(x) formatC(x, digits = digits, format = "f"))
  esc <- function(x) gsub("_", "\\_", as.character(x), fixed = TRUE)
  header <- paste(esc(names(df)), collapse = " & ")
  body <- apply(df, 1L, function(x) paste(esc(x), collapse = " & "))
  alignment <- paste(rep("l", ncol(df)), collapse = "")
  if (nrow(df) > 35L) {
    lines <- c(
      "\\begingroup", "\\scriptsize",
      paste0("\\begin{longtable}{", alignment, "}"),
      paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
      "\\toprule", paste0(header, " \\\\"), "\\midrule", "\\endfirsthead",
      "\\toprule", paste0(header, " \\\\"), "\\midrule", "\\endhead",
      paste0(body, " \\\\"),
      "\\bottomrule", "\\end{longtable}", "\\endgroup"
    )
  } else {
    lines <- c(
      "\\begin{table}[!htbp]", "\\centering", "\\small",
      paste0("\\caption{", caption, "}"),
      paste0("\\label{", label, "}"),
      paste0("\\begin{tabular}{", alignment, "}"),
      "\\toprule", paste0(header, " \\\\"), "\\midrule",
      paste0(body, " \\\\"),
      "\\bottomrule", "\\end{tabular}", "\\end{table}"
    )
  }
  writeLines(lines, path)
}

effect_result <- NULL
geometry_result <- NULL
if (study %in% c("both", "effect")) effect_result <- run_effect_study()
if (study %in% c("both", "geometry")) geometry_result <- run_geometry_study()

if (!is.null(effect_result)) {
  saveRDS(effect_result, file.path(generated_dir, paste0(run_id, "_effect_modification.rds")))
  if (!isTRUE(effect_result$pilot_only)) {
    write.csv(effect_result$h_summary,
              file.path(table_dir, "effect_modification_h_summary.csv"), row.names = FALSE)
    write.csv(effect_result$causal_summary,
              file.path(table_dir, "effect_modification_causal_summary.csv"), row.names = FALSE)
    write.csv(effect_result$label_summary,
              file.path(table_dir, "effect_modification_label_summary.csv"), row.names = FALSE)
    p <- make_effect_figure(effect_result)
    save_ggplot_grid(
      p, file.path(figure_dir, "robustness_effect_modification.png"),
      widths = c(1), heights = c(0.85, 1),
      layout_matrix = matrix(c(1, 2), nrow = 2, byrow = TRUE),
      width = 12, height = 7
    )
    save_ggplot_grid(
      p, file.path(figure_dir, "robustness_effect_modification.pdf"),
      widths = c(1), heights = c(0.85, 1),
      layout_matrix = matrix(c(1, 2), nrow = 2, byrow = TRUE),
      width = 12, height = 7
    )
  }
}

if (!is.null(geometry_result)) {
  saveRDS(geometry_result, file.path(generated_dir, paste0(run_id, "_geometry_transport.rds")))
  if (!isTRUE(geometry_result$pilot_only)) {
    write.csv(geometry_result$causal_summary,
              file.path(table_dir, "geometry_transport_summary.csv"), row.names = FALSE)
    write.csv(geometry_result$label_summary,
              file.path(table_dir, "geometry_transport_label_summary.csv"), row.names = FALSE)
    p <- make_geometry_figure(geometry_result)
    save_ggplot_grid(
      p, file.path(figure_dir, "robustness_geometry_transport.png"),
      widths = c(1), heights = c(0.85, 1),
      layout_matrix = matrix(c(1, 2), nrow = 2, byrow = TRUE),
      width = 12, height = 7
    )
    save_ggplot_grid(
      p, file.path(figure_dir, "robustness_geometry_transport.pdf"),
      widths = c(1), heights = c(0.85, 1),
      layout_matrix = matrix(c(1, 2), nrow = 2, byrow = TRUE),
      width = 12, height = 7
    )
  }
}

config <- list(
  study = study, test_mode = test_mode, pilot_only = pilot_only,
  n_rep = n_rep, pilot_reps = pilot_reps, n_truth = n_truth, n_cf = n_cf,
  n_workers = n_workers, h_sem_rounds = h_rounds,
  base_seed = base_seed, observed_design_seed_effect = base_seed + 101L,
  path_seed_geometry = base_seed + 201L,
  observed_design_seed_geometry = base_seed + 202L,
  conditioning = "replication-specific truth and estimates conditional on the same pre-treatment history",
  effect_design = "32x3 realizations; one effect-aware SEM + naive/oracle heterogeneous fits per realization",
  geometry_design = "simulate, SEM, and fit under each coarseness-path allocation; evaluate fixed all-or-nothing DAITE",
  grid = c(10L, 10L), treatment_time = treatment_time,
  observation_window = c(0, end_time), control_params = control_params,
  treated_params = treated_params, sem_control = sem_control
)
saveRDS(
  list(config = config, effect = effect_result, geometry = geometry_result),
  file.path(generated_dir, paste0(run_id, ".rds"))
)
message("[structured] outputs: ", generated_dir)
