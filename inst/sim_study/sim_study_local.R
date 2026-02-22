# This script contains examples of applying the full "Kresins-Clark" estimation procedure (WHATEVER THIS PROCEDURE ACTUALLY IS ! :))
library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(pbapply)
library(parallel)
library(doParallel)
library(R.utils)
library(caret)

# =================================
# GLOBAL PARAMS:
SAVE_DIR = '/Users/dac6/Library/CloudStorage/GoogleDrive-clark.a.duncan@gmail.com/My Drive/my_laptop/Documents/Academics/computation/r_code/CausalPointProcess'
# Data Geeneration Process
TREAT_PROP = 0.5
OMEGA = c(0,10^2,0,10^2)
END_TIME = 110
TREATMENT_TIME <- 10
NX = 10
NY = 10

# RUNTIME PARAMS
N_CORES <- 7
SIM_SIZE = N_CORES # number of observed data to realize
N_SIMS = 10 # number of simulations to estimate the ATE
N_TAU_SIMS = 10 # number of simulations to estimate the tau_i with
N_TAU_I = 10 # number of tau_i to estimate the ATE
N_TAU_I_TRUE = 100
N_PROPOSALS <- 1000
EM_ITER = 1000
SEM_EM_ADAPTIVE_ITER <- 1000
MAX_TIME <- 10000*(END_TIME*OMEGA[2]*OMEGA[4]/(1e6))
TIME_INT <- END_TIME - TREATMENT_TIME

# =====================================
# Simulation Stduy Description
# 1) Hawkes - Hawkes
# =====================================

make_cluster <- function(N_CORES){
  # setup the cluster:
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  estimatedparams<-list()
  # export libraries to cluster:
  clusterEvalQ(cl, {
    library(data.table)
    library(dplyr)
    library(spatstat)
    library(cubature)
    library(R.utils)
    library(caret)
    source('R/hawkes_extra.R')
    source('R/utils.R')
    source('R/labelings.R')
    source('R/thinings.R')
    source('R/em_algorithm.R')
  })
  clusterExport(cl, c("OMEGA",
                      "partition",
                      "N_SIMS",
                      "MAX_TIME",
                      "N_TAU_I",
                      "N_TAU_SIMS",
                      "TREATMENT_TIME",
                      "END_TIME",
                      "MAX_TIME",
                      "TIME_INT",
                      "tilenames",
                      "state_spaces",
                      "partition_processes",
                      "treated_partitions",
                      "treated_state_space",
                      "control_state_space",
                      "hawkes_par_1",
                      "hawkes_par_2",
                      "N_PROPOSALS",
                      "EM_ITER",
                      "SEM_EM_ADAPTIVE_ITER"))
  return(cl)
}

hawkes_par_1 <- list(mu = 8,
                     alpha = 0.01,
                     beta = 10,
                     K = 0.8)
hawkes_par_2 <- hawkes_par_1
hawkes_par_2$K <- 0.2

# True all/nothing ATE
partition <- quadrats(X=OMEGA, nx = NX, ny = NY)
all_nothing_ATE <- (hawkes_par_2$mu*TIME_INT*(1/(1-hawkes_par_2$K)) -
                      hawkes_par_1$mu*TIME_INT*(1/(1-hawkes_par_1$K)))/partition$n
partition_processes <- rep("control",partition$n)
partition_processes[sample(1:(NX*NY),NX*NY*TREAT_PROP)] <- "treated"
treated_idx <- partition_processes == "treated"
control_state_space <- as.owin(partition[!treated_idx])
treated_state_space <- as.owin(partition[treated_idx])
treated_partitions = tilenames(partition[which(partition_processes == "treated")])
if(partition_processes[1]=="treated"){
  state_spaces <- list(treated_state_space,control_state_space)
}else{
  state_spaces <- list(control_state_space,treated_state_space)
}

# get the true average 1-flip:
cl = make_cluster(N_CORES)
t<-proc.time()
tau_i_estim <- parSapply(1:partition$n,cl = cl,FUN = function(i){
  result <- tau_i(i,
                  partition = partition,
                  treated_partitions = treated_partitions,
                  statespace = OMEGA,
                  windowT =  c(TREATMENT_TIME,END_TIME),
                  control_pp = hawkes_par_1,
                  treated_pp = hawkes_par_2,
                  n_sim = N_TAU_I_TRUE)
  result
})
t_true_one_flip <- proc.time()[3] - t[3]
print(t_true_one_flip)
print(summary(tau_i_estim))
true_tau_1 <- mean(tau_i_estim)

# Generate some observed data
t <- proc.time()
obs_data <- lapply(1:SIM_SIZE,function(i){
  pre_treat <- sim_hawkes(params = hawkes_par_1,
                          windowT = c(0,TREATMENT_TIME),
                          windowS = OMEGA,
                          background_realization = NULL,
                          filtration = NULL,
                          optimized = TRUE)
  pre_treat$process <- "control"
  pre_treat$location_process <- "control"

  combined <- generate_inhomogeneous_hawkes(
    Omega = OMEGA,
    partition = partition,
    time_window = c(TREATMENT_TIME,END_TIME),
    partition_processes = partition_processes,
    state_spaces = state_spaces,
    filtration = pre_treat,
    hawkes_params = list(control = hawkes_par_1,treated = hawkes_par_2)
  )
  return(as.data.frame(combined))
})
t_obs_data <- (proc.time() - t)[3]
print(paste("Generating ",SIM_SIZE," simulations took ",t_obs_data," seconds."))

# Check number of points vs expected
tmp <- sapply(obs_data,function(x){sum(x$process == 'treated')})
hist(tmp,xlab ="Observed number of points in treatment area",main = 'Treated')
abline(v = hawkes_par_2$mu*TIME_INT*mean(partition_processes =='treated')*(1/(1-hawkes_par_2$K)),col = "blue", lwd = 2, lty = 2)
abline(v = mean(tmp),col = "red", lwd = 2, lty = 2)
legend("topright",legend = c("Expected","mean"),col = c("blue","red"),lty = c(2,2),lwd= c(2,2))

tmp <- sapply(obs_data,function(x){sum(x$t < TREATMENT_TIME)})
hist(tmp,xlab ="Observed number of points in control area",main = 'Control')
expected <- hawkes_par_1$mu*TREATMENT_TIME*(1/(1-hawkes_par_1$K))
abline(v = expected,col = "blue", lwd = 2, lty = 2)
abline(v = mean(tmp),col = "red", lwd = 2, lty = 2)
legend("topright",legend = c("Expected","mean"),col = c("blue","red"),lty = c(2,2),lwd= c(2,2))

# Check looks reasonable
as.data.frame(obs_data[[1]]) %>% filter(t < TREATMENT_TIME) %>% group_by(process,location_process) %>% summarize(n())
as.data.frame(obs_data[[1]]) %>% filter(t >= TREATMENT_TIME) %>% group_by(process,location_process) %>% summarize(n())

# Plot the point pattern with different colors
plot_1 <- plot_pp(as.data.frame(obs_data[[1]]) %>% filter(t > TREATMENT_TIME),
        partition = partition,
        title = "Simulated Hawkes-Hawkes Process - Complete Evolution"
)
print(plot_1)
ggsave('simulated_hawkes_hawkes_process.pdf',plot_1)

plot_2 <- plot_pp(as.data.frame(obs_data[[1]]) %>% filter(t < TREATMENT_TIME),
        partition = partition,
        title = "Simulated Hawkes-Hawkes Process - Pre Treatment"
)
print(plot_2)
ggsave('simulated_hawkes_hawkes_process_pretreat.pdf',plot_2)


# ==========================
# Labelling
# If we knew the labels we would be fine
# but we don't - but we can guess
# ==========================

pp_labeled_oracle <- lapply(obs_data,function(s){
  pre <- as.data.frame(s) %>% filter(t < TREATMENT_TIME) %>% mutate(inferred_process = "control",location_process = "control")
  post <- as.data.frame(s) %>% filter(t >= TREATMENT_TIME)
  return(rbind(pre,oracle_labeling(post)))
})

pp_labeled_naive <- lapply(obs_data,function(s){
  pre <- as.data.frame(s) %>% filter(t < TREATMENT_TIME) %>% mutate(inferred_process = "control",location_process = "control")
  post <- as.data.frame(s) %>% filter(t >= TREATMENT_TIME)
  return(rbind(pre,naive_labeling(post)))
})

# Get some labelling proposals with a reasonable sampler.
# Don't bother trying to select one of these - since we can't ever basically get a good one
cl <- make_cluster(N_CORES)
t<- proc.time()
labelling_proposals <- parLapply(cl=cl,X=obs_data,fun = function(s){
  pre <- as.data.frame(s) %>% filter(t < TREATMENT_TIME)
  post <- as.data.frame(s) %>% filter(t >= TREATMENT_TIME)
  pre$location_process <- "control"
  pre$inferred_process <- NULL
  lapply(1:N_PROPOSALS,function(i){
    tmp <- simulation_labeling_hawkes_hawkes(post,
                                             partition = partition,
                                             partition_process = partition_processes,
                                             statespace = OMEGA,
                                             state_spaces = state_spaces,
                                             windowT = c(TREATMENT_TIME,END_TIME),
                                             hawkes_params_control = hawkes_par_1,
                                             hawkes_params_treated = NULL,
                                             change_factor = 1,
                                             filtration = pre
    )
    pre$inferred_process <- "control"
    rbind(pre,tmp)
  })
})
t_labelling_props <-(proc.time()-t)[3]
print(paste0("Proposing ",
             length(labelling_proposals[[1]]),
             " travesals for each observation took ",
             t_labelling_props," seconds"))

# select the very best one as the best proposal:
proposal_accuracies <- lapply(labelling_proposals,function(x){
  tmp <- sapply(x,function(y){
    keep <- which(y$t > TREATMENT_TIME)
    acc <- confusionMatrix(as.factor(y$inferred_process[keep]),
                           as.factor(y$process[keep]))$overall[['Accuracy']]
    return(acc)})
  return(tmp)
})

pp_labeled_best_proposal <- lapply(1:length(labelling_proposals),function(i){
  best_proposal <- which.max(proposal_accuracies[[i]])
  return(labelling_proposals[[i]][[best_proposal]])
})

# funky EM style thing Duncan came up with
t<-proc.time()[3]
EM_max_style <- parLapply(cl = cl,X=obs_data,function(x){
  total_points <- sum(x$location_process == 'treated' & x$t >= TREATMENT_TIME)
  mu_start <- total_points / (TIME_INT)
  params_init = list(mu = mu_start,alpha = 0.1,beta = TIME_INT/10,K=0.1)
  t <- proc.time()
  # takes around 2 seconds per iteration
  tmp <- em_style_labelling(pp_data = x,
                     partition = partition,
                     partition_processes = partition_processes,
                     statespace = OMEGA,
                     state_spaces = state_spaces,
                     time_window = c(TREATMENT_TIME,END_TIME),
                     treatment_time = TREATMENT_TIME,
                     hawkes_params_control = hawkes_par_1,
                     hawkes_params_treated = params_init,
                     param_update_cadence = 10,
                     proposal_update_cadence = 1,
                     update_starting_data = T,
                     include_starting_data = F,
                     metric = 'post_likelihood',
                     optim_method = 'max',
                     iter = EM_ITER,
                     n_props = 10,
                     change_factor = 0.02,
                     MCMC_style = F,
                     verbose = F)
   print((proc.time() - t)[3])
  #tmp
  plot(tmp$metrics,tmp$accuracies)
  return(tmp)
})
t_em_post <- proc.time()[3] - t

pp_labelled_em_post <- lapply(EM_max_style,function(x){
  return(x$labelling)
})

# adaptive EM:
labellings <- list(oracle = pp_labeled_oracle,
                   naive = pp_labeled_naive,
                   best = pp_labeled_best_proposal,
                   EM_style = pp_labelled_em_post
                   )

# ================================
# Check labelling results :
# ===============================
class_metrics <- mapply(labellings,names(labellings),FUN=function(x,n){
  tmp <- lapply(x,function(y){
    keep <- which(y$t > TREATMENT_TIME)
    conf <- confusionMatrix(as.factor(y$inferred_process[keep]),
                           as.factor(y$process[keep]))
    # get the likelihood:
    realiz = y %>% filter(t >= TREATMENT_TIME)
    include = which(realiz$inferred_process == 'control')
    if(length(include) == 0){
      return(-Inf)
    }
    control_lik <- loglik_hawk(params = unlist(hawkes_par_1),
                               realiz = realiz[include,],
                               dists = NULL,
                               windowT = c(TREATMENT_TIME,END_TIME),
                               windowS = OMEGA,
                               zero_background_region = treated_state_space,
                               density_approx = FALSE,
                               numeric_integral= FALSE)

    include = which(realiz$inferred_process == 'treated')
    if(length(include) == 0){
      return(-Inf)
    }
    treat_lik <- loglik_hawk(params = unlist(hawkes_par_2),
                             realiz = realiz[include,],
                             dists = NULL,
                             windowT = c(TREATMENT_TIME,END_TIME),
                             windowS = OMEGA,
                             zero_background_region = control_state_space,
                             density_approx = FALSE,
                             numeric_integral= FALSE)
    tmp <- control_lik + treat_lik
    names(tmp) <- c("post_lik")
    return(c(conf$overall,conf$byClass,tmp))
  })
  tmp <- as.data.frame(do.call(rbind,tmp))
  tmp$name <- n
  tmp
},SIMPLIFY = F)
class_metrics <- do.call(rbind,class_metrics)
class_metrics <- class_metrics %>%
  group_by(name) %>%
  summarize(mean(Accuracy),
            mean(post_lik),
            mean(`Pos Pred Value`),
            mean(`Neg Pred Value`)) %>%
  ungroup()
class_metrics

# =====================================================================
# Perturb oracle labelling and see how it affects the different metrics
# =====================================================================
# test likelihood of the control process on perturbing some of the labels in
# the treatment region
if(FALSE){
  k  <- 1
  max_pert <- 20
  samp <- 20
  dat <- pp_labeled_oracle[[k]]
  t <- proc.time()
  cl <- make_cluster(N_CORES)
  clusterExport(cl,c("max_pert","samp","dat","control_state_space"))
  oracle_perturb <- parLapply(cl=cl,X = 1:max_pert,fun = function(p){
    result <- lapply(1:samp,function(j){
      set.seed(j)
      perturbed <- dat
      changes <- sample(which(dat$t > TREATMENT_TIME),p)
      current_labels <- dat$inferred_process[changes]
      new_labels <- c("control","treated")[1*(current_labels == "control") + 1]
      perturbed$inferred_process[changes] <- new_labels

      # fit the params on control data:
      post <- perturbed[perturbed$t >= TREATMENT_TIME,]
      pre <- perturbed[perturbed$t < TREATMENT_TIME,]
      control_param <- fit_hawkes(params_init = unlist(list(mu = length(post$t)/((END_TIME-TREATMENT_TIME)*area(control_state_space)/area(OMEGA)),
                                                            alpha = 1/OMEGA[2],# make it so that points within the region are fair game
                                                            beta = 1/END_TIME,
                                                            K = 0.9)),
                                  realiz = post[post$inferred_process == 'control',],
                                  windowT = c(TREATMENT_TIME,END_TIME),
                                  windowS = OMEGA,
                                  trace = 0,
                                  maxit = 1000,
                                  density_approx = FALSE,
                                  numeric_integral= FALSE,
                                  poisson_flag = FALSE,
                                  zero_background_region = treated_state_space,
      )

      pre_lik_c <- loglik_hawk(params = unlist(control_param$par),
                           realiz = pre,
                           windowT = c(0,TREATMENT_TIME),
                           windowS = OMEGA,
                           zero_background_region = NULL,
                           density_approx = FALSE,
                           numeric_integral= FALSE)

      post_lik_c_est <- loglik_hawk(params = unlist(control_param$par),
                               realiz = post %>% filter(inferred_process == 'control'),
                               windowT = c(TREATMENT_TIME,END_TIME),
                               windowS = OMEGA,
                               zero_background_region = treated_state_space,
                               density_approx = FALSE,
                               numeric_integral= FALSE)

      post_lik_c_true <- loglik_hawk(params = unlist(hawkes_par_1),
                                    realiz = post %>% filter(inferred_process == 'control'),
                                    windowT = c(TREATMENT_TIME,END_TIME),
                                    windowS = OMEGA,
                                    zero_background_region = treated_state_space,
                                    density_approx = FALSE,
                                    numeric_integral= FALSE)

      post_lik_t_true <- loglik_hawk(params = unlist(hawkes_par_2),
                                    realiz = post %>% filter(inferred_process == 'treated'),
                                    windowT = c(TREATMENT_TIME,END_TIME),
                                    windowS = OMEGA,
                                    zero_background_region = control_state_space,
                                    density_approx = FALSE,
                                    numeric_integral= FALSE)

      return(list(pre_lik_c = pre_lik_c,
             post_lik_c_est = post_lik_c_est,
             post_lik_c_true = post_lik_c_true,
             post_lik_t_true = post_lik_t_true,
             points_changed = p))
    })
    return(result)
  })
  print(proc.time() - t)

  pre_lik_c <- unlist(lapply(oracle_perturb,function(x){sapply(x,function(x){x$pre_lik_c})}))
  post_lik_c_est <- unlist(lapply(oracle_perturb,function(x){sapply(x,function(x){x$post_lik_c_est})}))
  post_lik_c_true <- unlist(lapply(oracle_perturb,function(x){sapply(x,function(x){x$post_lik_c_true})}))
  post_lik_t_true <- unlist(lapply(oracle_perturb,function(x){sapply(x,function(x){x$post_lik_t_true})}))
  post_total_lik_true <- unlist(lapply(oracle_perturb,function(x){sapply(x,function(x){x$post_lik_c_true + x$post_lik_t_true})}))

  all_vals <- unlist(oracle_perturb)

  # Pre lik metric
  plot(unlist(lapply(1:max_pert,function(i){rep(i,samp)})),
       pre_lik_c,
       pch = 19,
       col = "blue",
       xlab = "Number of points changed",
       ylab = "Log Likelihood",
       ylim = quantile(pre_lik_c,c(0.01,1)),
  )
  lines(1:max_pert,apply(sapply(oracle_perturb,function(x){sapply(x,function(x){x$pre_lik_c})}),2,mean),col = "red")

  # Post like metric:
  plot(unlist(lapply(1:max_pert,function(i){rep(i,samp)})),
       post_lik_c_true,
       pch = 19,
       col = "blue",
       xlab = "Number of points changed",
       ylab = "Log Likelihood",
       ylim = quantile(post_lik_c_true,c(0.01,1)),
  )
  lines(1:max_pert,apply(sapply(oracle_perturb,function(x){sapply(x,function(x){x$post_lik_c_true})}),2,mean),col = "red")

  plot(unlist(lapply(1:max_pert,function(i){rep(i,samp)})),
       post_total_lik_true,
       pch = 19,
       col = "blue",
       xlab = "Number of points changed",
       ylab = "Log Likelihood",
       ylim = quantile(post_total_lik_true,c(0.01,1)),
  )
  lines(1:max_pert,apply(sapply(oracle_perturb,function(x){sapply(x,function(x){x$post_lik_c_true + x$post_lik_t_true})}),2,mean),col = "red")
}

# ====================================
# Run the adaptive EM procedure:
# ====================================

t_em_pure <- proc.time()[3]
EM_results <- parLapply(cl=cl,X=obs_data,fun = function(dat){
  total_points <- sum(dat$location_process == 'treated' & dat$t >= TREATMENT_TIME)
  mu_start <- total_points / (TIME_INT)
  params_init = list(mu = mu_start,alpha = 0.1,beta = TIME_INT/10,K=0.1)

  result <- adaptive_SEM(pp_data = dat,
              partition = partition,
              partition_processes = partition_processes,
              statespace = OMEGA,
              time_window = c(TREATMENT_TIME,END_TIME),
              treatment_time = TREATMENT_TIME,
              hawkes_params_control = hawkes_par_1,
              hawkes_params_treated = params_init,
              N_labellings = 100,
              N_iter = 10,
              verbose = T,

              adaptive_control = list(
                param_update_cadence = 10,
                proposal_update_cadence = 1,
                state_spaces = NULL,
                iter  = SEM_EM_ADAPTIVE_ITER,
                n_props = 10,
                change_factor = 0.02,
                include_starting_data = F,
                update_starting_data = T,
                verbose = T)
              )
  return(result)
})
t_em_pure <- proc.time()[3] - t_em_pure

save.image("sim_study_results.RData")

# tuning:
# result <- adaptive_SEM(pp_data = dat,
#                        partition = partition,
#                        partition_processes = partition_processes,
#                        statespace = OMEGA,
#                        time_window = c(TREATMENT_TIME,END_TIME),
#                        treatment_time = TREATMENT_TIME,
#                        hawkes_params_control = hawkes_par_1,
#                        hawkes_params_treated = params_init,
#                        N_labellings = 1000,
#                        N_iter = 10,
#                        verbose = T,
#
#                        adaptive_control = list(
#                          param_update_cadence = 50,
#                          proposal_update_cadence = 1,
#                          state_spaces = NULL,
#                          # ====== Key tuning params ===========
#                          iter  = 5000,
#                          n_props = 10,
#                          change_factor = 0.02,
#                          # ====================================
#                          include_starting_data = F,
#                          update_starting_data = T,
#                          verbose = T)
# )

# ====================================
# Estimate the ATE from the labellings:
# ====================================

# flatten the labelings for parallelism
keep <- names(labellings)
labellings <- labellings[keep]
labellings[["EM_full"]] <- lapply(EM_results,function(x){
  x$adaptive$adaptive_labelling
})

tasks <- list()
for (j in seq_along(labellings)) {
  labelling <- labellings[[j]]
  labelling_name <- names(labellings)[j]

  for (i in seq_along(labelling)) {
    tasks[[length(tasks) + 1]] <- list(
      j = j,
      i = i,
      x = labelling[[i]] %>% filter(t >= TREATMENT_TIME),
      labelling_name = labelling_name,
      hawkes_params = NULL
    )
  }
}
print(paste0("number of tasks to run is ",length(tasks)))

# need to add the EM_full params onto this:
# labelling not used - so just add on observed data with naive labelling"
for (i in seq_along(EM_results)) {
  tmp <- obs_data[[i]]
  tmp$inferred_process <- tmp$location_process

  tasks[[length(tasks) + 1]] <- list(
    j = length(labellings) + 1,
    i = i,
    x = tmp %>% filter(t >= TREATMENT_TIME),
    labelling_name = labelling_name,
    hawkes_params = list(control = EM_results[[i]]$hawkes_params_control,
                         treated = EM_results[[i]]$hawkes_params_treated)
  )
}

# Perform parallel processing on all tasks
t <- proc.time()
task_function <- function(task) {
  j <- task$j
  i <- task$i
  x <- task$x
  labeling_name <- task$labeling_name

  # Perform the computation with timeout and error handling
  result <- tryCatch({
    withTimeout(
      expr = ATE_estim_hawkes(
        statespace = OMEGA,
        partition = partition,
        observed_data = x,
        treated_partitions = tilenames(partition[which(partition_processes == "treated")]),
        hawkes_params = task$hawkes_params,
        n_sims = N_SIMS,
        n_tau_sims = N_TAU_SIMS,
        n_tau_i = N_TAU_I,
        windowT = c(TREATMENT_TIME, END_TIME),
        windowS = OMEGA,
        maxit = 1000,
        poisson_flags = list(control = FALSE, treated = FALSE)
      ),
      timeout = MAX_TIME,
      onTimeout = "error"  # Throws an error if timeout occurs
    )
  }, TimeoutException = function(ex) {
    warning(paste("Labeling:", labeling_name, "- Iteration:", i,
                  "timed out after", MAX_TIME, "seconds"))
    return(NULL)  # Return NULL or any default value
  }, error = function(e) {
    warning(paste("Labeling:", labeling_name, "- Iteration:", i,
                  "error:", e$message))
    return(e)  # Handle other errors if necessary
  })
  # Return the result along with j and i for reconstruction
  return(list(j = j, i = i, result = result))
}
cl <- make_cluster(N_CORES)
results_flat <- parLapply(cl, tasks,fun = task_function)
# Calculate and print elapsed time
t4 <- (proc.time() - t)[3]
print(paste("Elapsed time:", t4, "seconds"))
stopCluster(cl)

# un flatten results_h_h
results_h_h <- vector("list", length(labellings) + 1)  # +1 for EM_full
names(results_h_h) <- c(names(labellings),"EM_full")

# Populate results_h_h based on results_flat
for (task_result in results_flat) {
  j <- task_result$j
  i <- task_result$i
  result <- task_result$result

  # Initialize sublist if it doesn't exist
  if (is.null(results_h_h[[j]])) {
    results_h_h[[j]] <- vector("list", length(obs_data))
    print(length(results_h_h[[j]]))
  }
  # Assign the result to the appropriate position
  results_h_h[[j]][[i]] <- result
}

# drop values that timed out i.e. that are NULL
print("Count number of nulls")
sapply(results_h_h,function(x){sum(sapply(x,is.null))})
results_h_h_orig <- results_h_h
results_h_h <- lapply(results_h_h,function(x){x[!sapply(x,is.null)]})
sapply(results_h_h,function(x){sum(sapply(x,is.null))})
print("count number of errors")
sapply(results_h_h,function(x){sum((sapply(x,function(y){!is.null(y$message)})))})
results_h_h <- lapply(results_h_h,function(x){
  null <- sapply(x,function(y){!is.null(y$message)})
  x[!null]
}
)
results_h_h <- results_h_h[sapply(results_h_h,length) > 0]
sapply(results_h_h,function(x){sum((sapply(x,function(y){!is.null(y$message)})))})
results_flat_h_h <- do.call(c,results_h_h)

# get the mean estimates of the hawkes process params in each case:
control_params <- lapply(results_h_h,function(x){
  tmp <- lapply(x,function(y){
    as.data.frame(y$control_pp)
  })
  rbindlist(tmp)
})
treated_params <- lapply(results_h_h,function(x){
  tmp <- lapply(x,function(y){
    as.data.frame(y$treated_pp)
  })
  rbindlist(tmp)
})

# ===================================
# Parameter Estimate Results
# ===================================
for(i in 1:length(control_params)){
  print(paste("Mean control Hawkes Params for ",names(control_params)[i]))
  print(colMeans(control_params[[i]]))
}
for(i in 1:length(treated_params)){
  print(paste("Mean treated Hawkes Params for ",names(treated_params)[i]))
  print(colMeans(treated_params[[i]]))
}

# make a boxplot for control param estimation for different labellings:
control_params_df <- mapply(control_params,names(control_params),FUN = function(df,n){
  df$labelling <- n
  return(df)
},SIMPLIFY = FALSE)
control_params_df <- do.call(rbind,control_params_df) %>%
  reshape2::melt(id.vars = "labelling") %>%
  rename("labelling" = "variable",
         "param" = "variable",
         "value" = "value")

control_param_plots <- list()
for(n in unique(control_params_df$param)){
  tmp <- control_params_df %>%
    filter(param == n)
  plot <- ggplot(tmp) +
    geom_boxplot(aes(x = labelling, y = value)) +
    labs(
      title = paste0(n," estimates for Hawkes-Hawkes"),
      x = "Method",
      y = "Parameter Estimate",
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(face = "bold")
    )+
    # add hline:
    geom_hline(yintercept = hawkes_par_1[[n]],linetype = "dashed",color = "red")
  control_param_plots[[n]] <- plot
}
# combine the grobs:
plot_3 <- gridExtra::grid.arrange(grobs = control_param_plots,ncol = 2)
print(g)
ggsave('hawkes_hawkes_control_param_estimates.pdf',plot_3)


# make a boxplot for treated param estimation for different labellings:
treated_params_df <- mapply(treated_params,names(treated_params),FUN = function(df,n){
  df$labelling <- n
  return(df)
},SIMPLIFY = FALSE)
treated_params_df <- do.call(rbind,treated_params_df) %>%
  reshape2::melt(id.vars = "labelling") %>%
  rename("labelling" = "variable",
         "param" = "variable",
         "value" = "value")

treated_param_plots <- list()
for(n in unique(treated_params_df$param)){
  tmp <- treated_params_df %>%
    filter(param == n)
  plot <- ggplot(tmp) +
    geom_boxplot(aes(x = labelling, y = value)) +
    labs(
      title = paste0(n," estimates for Hawkes-Hawkes example"),
      x = "Method",
      y = "Parameter Estimate",
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(face = "bold")
    )+
    # add hline:
    geom_hline(yintercept = hawkes_par_2[[n]],linetype = "dashed",color = "red")
  treated_param_plots[[n]] <- plot
}
# combine the grobs:
plot_4 <- gridExtra::grid.arrange(grobs = treated_param_plots,ncol = 2)
print(g)
ggsave('hawkes_hawkes_treated_param_estimates.pdf',plot_4)

# ===================================
# Estimand Results
# ===================================

# Get the estimates for each of the quantities from each of the procedures:
data_list <- mapply(results_flat,
                    tasks[1:length(results_flat)],
                    1:length(results_flat),
                    FUN = function(result,task,i){
                      if(is.null(result$result) | !is.null(result$result$message)){return(NULL)}
                      df = as.data.frame(result$result[c("all_nothing_sim",
                                                         "all_nothing_theory",
                                                         "tau_1_estim",
                                                         "ATE_total",
                                                         "ATE_treatment",
                                                         "ATE_spillover",
                                                         "ATE_naive")]
                      )

                      control_pp <- result$result$control_pp
                      treated_pp <- result$result$treated_pp
                      windowT <- c(TREATMENT_TIME, END_TIME)

                      df = data.frame(all_nothing_sim = mean(result$result$all_nothing_sim$ATE),
                                      all_nothing_theory = mean(result$result$all_nothing_theory$ATE),
                                      points_per_tile_control_theory = mean(result$result$all_nothing_theory$c_mean),
                                      points_per_tile_treated_theory = mean(result$result$all_nothing_theory$t_mean),
                                      tau_1_estim = result$result$tau_1_estim,
                                      ATE_total = result$result$ATE_total,
                                      ATE_treatment = result$result$ATE_treatment,
                                      ATE_spillover = result$result$ATE_spillover,
                                      ATE_naive = result$result$ATE_naive,
                                      labelling = task$labelling_name
                      )
                      return(df)
                    },SIMPLIFY = FALSE)
data_list <- data_list[!sapply(data_list,is.null)]
data <- do.call(rbind,data_list)
data <- reshape2::melt(data,id.vars = "labelling") %>%
  rename("labelling" = "variable",
         "method" = "variable",
         "ATE_estim" = "value")

# order the data
data$method <- factor(data$method, levels = c("all_nothing_sim",
                                              "all_nothing_theory",
                                              "points_per_tile_control_theory",
                                              "points_per_tile_treated_theory",
                                              "tau_1_estim",
                                              'ATE_total',
                                              'ATE_treatment',
                                              'ATE_spillover',
                                              'ATE_naive'))
data$estimand <- c("all-nothing",
                   "all-nothing",
                   "control_points",
                   "treated_points",
                   "one-flip-spillover",
                   "non-param",
                   "non-param",
                   "non-param",
                   "non-param"
)[data$method]

means <- data %>%
  group_by(method,labelling) %>%
  summarise(mean_ATE_estim = mean(ATE_estim))
data <- data %>% left_join(means,by = c("method","labelling"))
# Separate data frame for lines
lines_data <- data.frame(
  method = data$method,
  labelling = data$labelling,
  mean_ATE_estim = data$mean_ATE_estim,
  true_1_flip = true_tau_1,
  all_nothing_ATE = unique(all_nothing_ATE)  # Replace with the actual all_nothing_ATE value if it's the same for all
) %>% distinct()

data$labelling<- factor(data$labelling,
                        # get the ordering the same as before
                        levels = c("oracle", "naive", "best", "EM_style", "EM_full")
                        )

# ====================
# All nothing boxplot
# ====================
# make a boxplot for the all_nothing estimate
# facet my simulation vs theoretical

tmp <- data %>% filter(method == "all_nothing_sim" | method == "all_nothing_theory")
tmp$type <- sapply(levels(tmp$method)[tmp$method],function(x){if(grepl("sim",x)){"sim"}else{"theory"}})

# make boxplot
plot_5 <- ggplot(tmp %>% filter(type == 'theory')) +
  geom_boxplot(aes(x = labelling, y = ATE_estim)) +
  geom_hline(data = lines_data, aes(yintercept = all_nothing_ATE, linetype = "True All/Nothing ATE"),col = scales::hue_pal()(3)[1],linewidth = 1) +
  scale_linetype_manual(name = "", values = c("True All/Nothing ATE" = "solid",
                                              "Average Single Flip ATE" = "solid")) +
  labs(
    title = "All-Nothing ATE estimates for Hawkes-Hawkes example",
    x = "Method",
    y = "All-Nothing ATE Estimate",
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(all_nothing_ATE*1.5,3)+
  theme(legend.position = "bottom")
print(plot)

ggsave("Hawkes_Hawkes_all_nothing_ATE.pdf",plot_5)

# make boxplot for the all-nothing situation showing control and treatment effects:
tmp <- data %>% filter(method == "points_per_tile_control_theory" | method == "points_per_tile_treated_theory")
tmp$method <- droplevels(tmp$method)
tmp$type <- sapply(levels(tmp$method)[tmp$method],function(x){if(grepl("sim",x)){"sim"}else{"theory"}})

oracle_means <- tmp %>%
  filter(type == 'theory', labelling == "oracle") %>%
  group_by(method) %>%
  summarise(mean_ATE = median(ATE_estim))
true_means <- data.frame(method = c("points_per_tile_control_theory","points_per_tile_treated_theory"),
                         mean_ATE = c(hawkes_par_1$mu*TIME_INT*(1/(1-hawkes_par_1$K))/partition$n,
                                      hawkes_par_2$mu*TIME_INT*(1/(1-hawkes_par_2$K))/partition$n)
                         )
# make boxplot
plot_6 <- ggplot(tmp %>% filter(type == "theory"),
                 aes(x = labelling, y = ATE_estim, fill = method)) +
  geom_boxplot() +
  facet_wrap(~ method) +

  # map Oracle lines to colour = "Oracle mean"
  geom_hline(aes(yintercept = mean_ATE, colour = "Oracle mean"),
             data     = oracle_means,
             linetype = "dashed",
             size     = 1) +
  # map True   lines to colour = "True mean"
  geom_hline(aes(yintercept = mean_ATE, colour = "True mean"),
             data     = true_means,
             linetype = "dashed",
             size     = 1) +

  # now supply the colours and labels for those two levels
  scale_colour_manual(
    name   = "Reference",
    values = c(
      "Oracle mean" = "#0072B2",
      "True mean"   = "#D55E00"
    )
  ) +

  # keep your box-fill legend separate
  guides(
    fill   = guide_legend(order = 1),
    colour = guide_legend(order = 2)
  ) +

  labs(
    title = "Points-Per-Tile Hawkes-Hawkes example",
    x     = "Method",
    y     = "Points Per Tile"
  ) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  ylim(0, max(true_means$mean_ATE))

plot_6

print(plot_6)
# save the plot:
ggsave("Hawkes_Hawkes_points_per_tile.pdf",plot_6)


# ====================
# One flip boxplot
# ====================

tmp <- data %>% filter(method == "tau_1_estim" | method == "ATE_naive")
plot_7 <- ggplot(tmp) +
  geom_boxplot(aes(x = labelling, y = ATE_estim)) +
  geom_hline(data = lines_data, aes(yintercept = true_1_flip , linetype = "Average Single Flip ATE"),col = scales::hue_pal()(3)[1],linewidth = 1) +
  scale_linetype_manual(name = "", values = c("All/Nothing ATE" = "solid",
                                              "Average Single Flip ATE" = "solid")) +
  labs(
    title = "Average Single Flip ATE for Hawkes-Hawkes example",
    x = "Method",
    y = "Single Flip ATE Estimate",
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(min(tmp$ATE_estim)*0.9,all_nothing_ATE*1.2)+
  theme(legend.position = "bottom")
print(plot_7)

# save the plot:
ggsave("Hawkes_Hawkes_one_flip.pdf",plot_7)


for(i in 1:length(results_h_h$oracle)){
  print(i)
  # print the expected results
  treated_pp <- results_h_h$oracle[[i]]$treated_pp
  control_pp <- results_h_h$oracle[[i]]$control_pp
  expected_control <- control_pp$mu*(TIME_INT)*(1/(1-control_pp$K))/partition$n
  expected_treated <- treated_pp$mu*(TIME_INT)*(1/(1-treated_pp$K))/partition$n
  expected_ATE <- expected_treated - expected_control
  print(paste("Expected ATE for oracle sim ",i,": ",round(expected_ATE,2)))
  # print the simulated results:
  print(paste("Simulated ATE for oracle sim ",i,": ",round(mean(results_h_h$oracle[[i]]$ATE_total),2)))
}

# print what took the most time:
print(data.frame('true one flip calc' = t_true_one_flip,
                 'sim obs data'= t_obs_data,
                 'proposals' = t_labelling_props,
                 'em_style' = t_em_post,
                 'em_pure' = t_em_pure,
                 'ATE_estim_hawkes' = t4))

# Need to save the relevant things!
save.image("sim_study_results.RData")




