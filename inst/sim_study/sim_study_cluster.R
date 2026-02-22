# This simulation study is similar to the examples
# - reduce number of tiles
# - superposition of two poisson processes
# - hawkes-hawkes deviance from the oracle thinning

# This is the script that is run on cluster, sim_study_process.R is
# the script that actually makes the paper outputs.

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

# GLOBAL PARAMS:
TREAT_PROP = 0.5
N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 16))
SIM_SIZE = 100
N_SIMS = 10 # number of simulations to estimate the ATE
N_TAU_SIMS = 10 # number of simulations to estimate the tau_i with
N_TAU_I = 10 # number of tau_i to estimate the ATE
OMEGA = c(0,10^2,0,10^2)
END_TIME = 110
TREATMENT_TIME <- 10
MAX_TIME <- 10000*(END_TIME*OMEGA[2]*OMEGA[4]/(1e6))
N_PROPOSALS <- 1000
EM_ITER = 100
SAVE_DIR = '/sim_study_outputs'
NX = 10
NY = 10

TIME_INT <- END_TIME - TREATMENT_TIME

# =====================================
# Sim Study Description
# 2) Inhomogeneous Poisson with spatial dispersion kernel
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
    source('hawkes_extra.R')
    source('utils.R')
    source('labelings.R')
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
                      "EM_ITER"))
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
                  n_sim = 100)
  result
})
print(proc.time()-t)
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
print(paste("Generating ",SIM_SIZE," simulations took ",(proc.time()-t)[3]," seconds."))




