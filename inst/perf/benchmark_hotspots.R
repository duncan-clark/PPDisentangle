#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(spatstat.geom)
  library(spatstat.random)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) < 1L) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}
has_flag <- function(name) name %in% args

repo_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1], fixed = FALSE), mustWork = FALSE)), "..", ".."), mustWork = FALSE)
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)

quick <- has_flag("--quick")
run_decode <- has_flag("--run-decode")
out_file <- arg_value("--out", file.path(repo_root, "inst", "perf", "benchmark_results.json"))
label <- arg_value("--label", "run")
seed <- as.integer(arg_value("--seed", "123"))
if (!is.finite(seed)) seed <- 123L
sem_iter_arg <- as.integer(arg_value("--sem-iter", NA_character_))
sem_props_arg <- as.integer(arg_value("--sem-props", NA_character_))
set.seed(seed)

elapsed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}

mk_test_data <- function(seed_local = 42L) {
  set.seed(seed_local)
  omega <- c(0, 50, 0, 50)
  end_time <- 55
  treat_time <- 5
  nx <- 5
  ny <- 5
  hawkes_par_1 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.5)
  hawkes_par_2 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.2)
  partition <- spatstat.geom::quadrats(X = omega, nx = nx, ny = ny)
  partition_processes <- rep("control", partition$n)
  partition_processes[sample(seq_len(partition$n), floor(partition$n * 0.5))] <- "treated"
  treated_idx <- partition_processes == "treated"
  control_ss <- spatstat.geom::as.owin(partition[!treated_idx])
  treated_ss <- spatstat.geom::as.owin(partition[treated_idx])
  state_spaces <- if (partition_processes[1] == "treated") list(treated_ss, control_ss) else list(control_ss, treated_ss)

  pre <- sim_hawkes(hawkes_par_1, c(0, treat_time), omega)
  pre_df <- data.frame(
    x = pre$x, y = pre$y, t = pre$t,
    process = "control", location_process = "control",
    inferred_process = "control", background = pre$background
  )
  post <- generate_inhomogeneous_hawkes(
    omega, partition, c(treat_time, end_time),
    partition_processes, list(control = hawkes_par_1, treated = hawkes_par_2),
    state_spaces = state_spaces, filtration = pre_df
  )
  post_df <- as.data.frame(post)
  post_df$inferred_process <- post_df$location_process
  cols <- c("x", "y", "t", "process", "location_process", "inferred_process", "background")
  pre_df <- pre_df[, cols, drop = FALSE]
  post_df <- post_df[, cols, drop = FALSE]
  obs <- rbind(pre_df, post_df)
  obs$W <- 1

  list(
    omega = omega,
    end_time = end_time,
    treat_time = treat_time,
    partition = partition,
    partition_processes = partition_processes,
    state_spaces = state_spaces,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    obs = obs
  )
}

bench_likelihood <- function() {
  params <- c(mu = 30, alpha = 0.5, beta = 5, K = 0.3)
  window_t <- c(0, 10)
  window_s <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  zero_region <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  sim <- sim_hawkes(as.list(params), window_t, window_s)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))
  n_rep <- if (quick) 20L else 80L

  t_no_precomp <- elapsed(for (i in seq_len(n_rep)) {
    loglik_hawk_fast(params, realiz, window_t, window_s, zero_background_region = zero_region)
  })
  pc <- precompute_loglik_args(realiz, window_s, zero_region)
  t_precomp <- elapsed(for (i in seq_len(n_rep)) {
    loglik_hawk_fast(params, realiz, window_t, window_s,
      precomp = list(active_area = pc$active_area, in_zero_bg = pc$in_zero_bg_all))
  })

  ll_no <- loglik_hawk_fast(params, realiz, window_t, window_s, zero_background_region = zero_region)
  ll_pre <- loglik_hawk_fast(params, realiz, window_t, window_s,
    precomp = list(active_area = pc$active_area, in_zero_bg = pc$in_zero_bg_all))

  list(
    name = "hawkes_loglik_precomp",
    primary_time_sec = unname(t_precomp),
    no_precomp_sec = unname(t_no_precomp),
    precomp_sec = unname(t_precomp),
    internal_speedup = as.numeric(t_no_precomp / t_precomp),
    ll_abs_diff = as.numeric(abs(ll_no - ll_pre)),
    correctness_pass = is.finite(ll_no) && is.finite(ll_pre) &&
      (abs(ll_no - ll_pre) <= 20 || abs(ll_no - ll_pre) / max(1, abs(ll_no)) <= 0.01)
  )
}

bench_simulation <- function() {
  params <- list(mu = 40, alpha = 0.2, beta = 8, K = 0.35)
  window_t <- c(0, 8)
  window_s <- spatstat.geom::owin(xrange = c(0, 15), yrange = c(0, 15))
  n_rep <- if (quick) 30L else 120L

  set.seed(seed + 1L)
  t_fast <- elapsed({
    n_fast <- numeric(n_rep)
    for (i in seq_len(n_rep)) {
      s <- sim_hawkes_fast(params, window_t, window_s)
      n_fast[i] <- length(s$t)
    }
    assign(".n_fast", n_fast, envir = .GlobalEnv)
  })
  n_fast <- get(".n_fast", envir = .GlobalEnv); rm(".n_fast", envir = .GlobalEnv)

  set.seed(seed + 1L)
  t_ref <- elapsed({
    n_ref <- numeric(n_rep)
    for (i in seq_len(n_rep)) {
      s <- sim_hawkes(params, window_t, window_s)
      n_ref[i] <- length(s$t)
    }
    assign(".n_ref", n_ref, envir = .GlobalEnv)
  })
  n_ref <- get(".n_ref", envir = .GlobalEnv); rm(".n_ref", envir = .GlobalEnv)

  mu_diff <- abs(mean(n_fast) - mean(n_ref))
  sd_diff <- abs(sd(n_fast) - sd(n_ref))
  mu_rel <- mu_diff / max(1, abs(mean(n_ref)))
  sd_rel <- sd_diff / max(1, abs(sd(n_ref)))
  list(
    name = "hawkes_sim_fast_vs_ref",
    primary_time_sec = unname(t_fast),
    fast_sec = unname(t_fast),
    reference_sec = unname(t_ref),
    internal_speedup = as.numeric(t_ref / t_fast),
    mean_count_fast = as.numeric(mean(n_fast)),
    mean_count_ref = as.numeric(mean(n_ref)),
    sd_count_fast = as.numeric(sd(n_fast)),
    sd_count_ref = as.numeric(sd(n_ref)),
    mean_count_abs_diff = as.numeric(mu_diff),
    sd_count_abs_diff = as.numeric(sd_diff),
    mean_count_rel_diff = as.numeric(mu_rel),
    sd_count_rel_diff = as.numeric(sd_rel),
    correctness_pass = is.finite(mu_rel) && is.finite(sd_rel) &&
      mu_rel <= 0.05 && sd_rel <= 0.35
  )
}

bench_sem <- function() {
  d <- mk_test_data(seed + 2L)
  n_iter <- if (is.finite(sem_iter_arg)) sem_iter_arg else if (quick) 8L else 25L
  n_props <- if (is.finite(sem_props_arg)) sem_props_arg else if (quick) 8L else 20L
  set.seed(seed + 3L)
  sem_elapsed <- elapsed({
    sem_res <- em_style_labelling(
      pp_data = d$obs,
      partition = d$partition,
      partition_processes = d$partition_processes,
      statespace = d$omega,
      time_window = c(d$treat_time, d$end_time),
      treatment_time = d$treat_time,
      hawkes_params_control = d$hawkes_par_1,
      hawkes_params_treated = d$hawkes_par_2,
      update_control_params = TRUE,
      param_update_cadence = 5,
      proposal_update_cadence = 1,
      iter = n_iter,
      n_props = n_props,
      change_factor = 0.02,
      state_spaces = d$state_spaces,
      verbose = FALSE
    )
    assign(".sem_res", sem_res, envir = .GlobalEnv)
  })
  sem_res <- get(".sem_res", envir = .GlobalEnv); rm(".sem_res", envir = .GlobalEnv)
  acc <- sem_res$accuracies
  list(
    name = "sem_inner_loop",
    primary_time_sec = unname(sem_elapsed),
    elapsed_sec = unname(sem_elapsed),
    n_iter = length(acc),
    final_accuracy = as.numeric(tail(acc, 1)),
    mean_accuracy = as.numeric(mean(acc)),
    all_metrics_finite = all(is.finite(sem_res$metrics)),
    correctness_pass = all(is.finite(acc)) &&
      all(acc >= 0 & acc <= 1) &&
      nrow(sem_res$labelling) == nrow(d$obs) &&
      all(is.finite(sem_res$metrics))
  )
}

bench_decode <- function() {
  env <- Sys.getenv()
  old <- c(
    OK_RUN_DECODE = env[["OK_RUN_DECODE"]],
    OK_RUN_SENSITIVITY = env[["OK_RUN_SENSITIVITY"]],
    OK_REPORT_FORMATS = env[["OK_REPORT_FORMATS"]],
    OK_SEM_INNER_ITER = env[["OK_SEM_INNER_ITER"]],
    OK_DECODE_ITER = env[["OK_DECODE_ITER"]],
    OK_MEMORY_SAFE = env[["OK_MEMORY_SAFE"]],
    OK_CORES = env[["OK_CORES"]]
  )
  on.exit({
    for (nm in names(old)) {
      if (!nzchar(old[[nm]])) Sys.unsetenv(nm) else Sys.setenv(structure(old[[nm]], names = nm))
    }
  }, add = TRUE)

  Sys.setenv(
    OK_RUN_DECODE = "true",
    OK_RUN_SENSITIVITY = "false",
    OK_REPORT_FORMATS = "",
    OK_SEM_INNER_ITER = if (quick) "5" else "20",
    OK_DECODE_ITER = if (quick) "15" else "60",
    OK_MEMORY_SAFE = "true",
    OK_CORES = "1"
  )

  cmd <- file.path(repo_root, "inst", "oklahoma", "oklahoma_analysis.R")
  decode_elapsed <- elapsed({
    code <- system2("Rscript", c(shQuote(cmd), "--quick-check"), stdout = TRUE, stderr = TRUE)
    assign(".decode_code", attr(code, "status") %||% 0L, envir = .GlobalEnv)
  })
  status <- get(".decode_code", envir = .GlobalEnv); rm(".decode_code", envir = .GlobalEnv)
  out_rds <- file.path(repo_root, "output", "oklahoma_results.rds")
  has_results <- file.exists(out_rds)
  decode_stage <- NA_real_
  if (has_results) {
    rr <- tryCatch(readRDS(out_rds), error = function(e) NULL)
    if (!is.null(rr) && !is.null(rr$timing_df)) {
      td <- rr$timing_df
      keep <- td$stage %in% c("fit_G_decode_bivariate", "fit_H_decode_bivariate_kde")
      if (any(keep)) decode_stage <- sum(td$elapsed_sec[keep], na.rm = TRUE)
    }
  }
  list(
    name = "decode_end_to_end_quickcheck",
    primary_time_sec = unname(decode_elapsed),
    script_elapsed_sec = unname(decode_elapsed),
    decode_stage_sec = as.numeric(decode_stage),
    exit_status = as.integer(status),
    correctness_pass = identical(status, 0L) && has_results
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

results <- list(
  meta = list(
    label = label,
    seed = seed,
    quick = quick,
    run_decode = run_decode,
    timestamp = as.character(Sys.time()),
    repo_root = repo_root
  ),
  components = list(
    likelihood = bench_likelihood(),
    simulation = bench_simulation(),
    sem = bench_sem()
  )
)
if (isTRUE(run_decode)) {
  results$components$decode <- bench_decode()
}

ok <- all(vapply(results$components, function(x) isTRUE(x$correctness_pass), logical(1)))
results$overall_correctness_pass <- ok

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE), out_file)
cat("Wrote benchmark file:", out_file, "\n")
if (!ok) {
  cat("One or more correctness checks failed.\n")
  quit(status = 2L)
}
quit(status = 0L)
