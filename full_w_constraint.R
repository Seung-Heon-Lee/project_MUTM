# prerequistite
library(igraph) # to make graph
library(mipmatch) # use distmat function to construct distance matrix
library(Rcplex) # to solve IP
library(dplyr) # for data wrangling
library(tidyverse)
library("lattice") # for plotting
# define function
full_w_constraint <- function (X_mat, dist_mat, t_ind, max_control = NULL, max_treat = NULL, matched_control = NULL, mom_covs = NULL, 
                               mom_tols = NULL, w_constraints = FALSE) {
  # Basic parameters
  number_of_treated_nodes <- sum(t_ind)
  number_of_control_nodes <- sum(!t_ind)
  total <- number_of_control_nodes + number_of_treated_nodes
  if(is.null(max_control)) {max_control <- number_of_control_nodes}
  if(is.null(max_treat)) {max_treat <- number_of_treated_nodes}
  if(is.null(matched_control)) {matched_control <- number_of_control_nodes}
  
  # Define nodes
  nodes <- c("source", 
             paste0("t", 1:number_of_treated_nodes), 
             paste0("c", 1:number_of_control_nodes), 
             "sink", 
             "overflow")
  
  from = c(rep("source", number_of_treated_nodes), 
           rep(paste0("t", 1:number_of_treated_nodes), each = number_of_control_nodes), # all treated nodes to control nodes
           paste0("t", 1:number_of_treated_nodes), # all treated nodes to overflow node
           rep(paste0("c", 1:number_of_control_nodes), 2))
  
  to = c(paste0("t", 1:number_of_treated_nodes), 
         rep(paste0("c", 1:number_of_control_nodes), number_of_treated_nodes),
         rep("overflow", number_of_treated_nodes), 
         rep("sink", number_of_control_nodes),
         rep("overflow", number_of_control_nodes))
  
  
  cost = c(rep(0,number_of_treated_nodes), # source to each treated
           as.vector(t(dist_mat)), # distance matrix
           rep(0,number_of_treated_nodes), # each treated to overflow
           rep(0,number_of_control_nodes), # each control to sink
           rep(0,number_of_control_nodes)) # each control to overflow
  
  capacity = c(rep(max_control, number_of_treated_nodes), # source to each treated
               rep(1, number_of_control_nodes * number_of_treated_nodes), # each treated to control
               rep(max_control-1, number_of_treated_nodes), # each treated to overflow
               rep(1, number_of_control_nodes), # each control to sink
               rep(max_treat - 1, number_of_control_nodes) # each control to overflow
  )
  
  # Define edges
  edges <- data.frame(from, to,
                      capacity,
                      cost)
  
  # Create the graph
  g <- graph_from_data_frame(edges, directed=TRUE)
  
  # Define the supply for each node
  sply_source <- number_of_treated_nodes*max_control
  sply_over <- -(sply_source - matched_control)
  sply_sink <- -matched_control
  V(g)$supply <- c(sply_source, rep(0, total), sply_over, sply_sink) # supply of each control and treated node would be variable
  names(V(g)$supply) <- V(g)$name
  
  # let's solve MIP
  # Coefficeints of the objective function,  variables are flow of each edge
  cvec <- cost
  # Construct the constraint matrix
  if(w_constraints == TRUE) {Amat <- matrix(0, nrow=length(V(g)) + 2* ncol(mom_covs), ncol=length(E(g)))}
  if(w_constraints == FALSE) {Amat <- matrix(0, nrow=length(V(g)), ncol=length(E(g)))}
  
  # for each node outgoing edge is 1, incoming edge is -1. If there's no edge than 0
  for(i in 1:length(V(g))) {
    node <- V(g)$name[i]
    outgoing_edges <- incident(g, node, mode="out")
    incoming_edges <- incident(g, node, mode="in")
    
    # For outgoing edges of ith node, subtract flow from node supply
    Amat[i, outgoing_edges] <- 1 
    
    # For incoming edges, add flow to node supply
    Amat[i, incoming_edges] <- -1
  }
  # w/ CB constarint for  2* length(ncol(mom_covs)) rows
  if(w_constraints == TRUE) {
    for(n in 1:ncol(mom_covs)) {
      for(i in 1:number_of_treated_nodes) {
        node <- V(g)$name[i+1]
        outgoing_edges <- incident(g, node, mode="out")
        cov <- mom_covs[,n]
        coef <- numeric()
        for (j in (number_of_treated_nodes + 1):length(cov)) {
          coef <- c(coef, cov[i] - cov[j])
        }
        coef <- c(coef, 0)
        Amat[length(V(g)) + 2*n-1, outgoing_edges] <- coef
        Amat[length(V(g)) + 2*n, outgoing_edges] <- -coef
      }
    }
  }
  # Vector of RHS of constraints (i.e, supply for each node, 0 if not specified)
  bvec <- V(g)$supply
  if(w_constraints == TRUE) {bvec <- c(bvec, rep(mom_tols, each = 2))}
  
  # Bounds of variables
  lb <- rep(0, length(E(g)))
  ub <- c(as.vector(capacity))
  
  # the optimization direction
  objsense = c('min')
  
  # The direction of the inequality in each constraint. 'E' stands for equal
  # 'E' for supply constraints, 'L' for CB constraints
  if(w_constraints == TRUE ) {sense = c(rep('E',length(V(g))), rep("L", 2 * ncol(mom_covs)))}
  if(w_constraints == FALSE ) {sense = c(rep('E',length(V(g))))}
  
  # Types of variables, Integer
  vtype <- "I"
  
  # Solve the model
  ptm = proc.time()
  result <- Rcplex(cvec = cvec, Amat = Amat, bvec = bvec, Qmat = NULL,
                   lb = lb, ub = ub, control = list(),
                   objsense = objsense, sense = sense, vtype = vtype, n = 1)
  time = (proc.time() - ptm)[3]
  
  # get solution of full matching problem
  edge_ids <- as_ids(E(g))
  optimized_flows <- result$xopt
  df <- data.frame(edges = edge_ids, optimized_flows = optimized_flows)
  df_filtered <- df %>%
    filter(grepl("^t\\d+\\|c\\d+$", edges) & df$optimized_flows == 1) %>%
    separate(edges, into = c("treated", "control"), sep = "\\|")
  
  # Group by treated node and combine connected control nodes into a list
  df_grouped <- df_filtered %>%
    group_by(treated) %>%
    summarise(control = list(control)) %>% 
    group_by(control) %>% 
    summarise(treated = list(treated)) %>%
    mutate(stratum = map2(control, treated, ~ paste(c(.x, .y), collapse = " ")))
  strata_mip <- lapply(df_grouped$stratum, function(x) unlist(strsplit(x, split = " ")))
  num_of_stratum <- length(strata_mip)
  
  # calculate differences in each stratum
  rownames(X_mat) <- c(paste0("t", 1:number_of_treated_nodes),paste0("c", 1:number_of_control_nodes))
  lapply(strata_mip, function(x) X_mat[x,])
  t_avg <- matrix(0, nrow = num_of_stratum, ncol = ncol(X_mat) + 1)
  c_avg <- matrix(0, nrow = num_of_stratum, ncol = ncol(X_mat) + 1)
  
  # calculate smd after matching
  for(i in 1:num_of_stratum) {
    size <- length(strata_mip[[i]])
    df <- X_mat[strata_mip[[i]],]
    matrix(df[startsWith(rownames(df), "t"), ], ncol = ncol(X_mat))
    t_avg[i,] <- c(colMeans(matrix(df[startsWith(rownames(df), "c"), ], ncol = ncol(X_mat))), size)
    c_avg[i,] <- c(colMeans(matrix(df[startsWith(rownames(df), "t"), ], ncol = ncol(X_mat))), size)
  }
  mean_diff <- matrix(0, nrow = ncol(X_mat), ncol = 2)
  for(i in 1:ncol(X_mat)) {
    for(j in 1:num_of_stratum) {
      mean_diff[i,1] <- mean_diff[i,1] + (t_avg[j,ncol(X_mat)+1] * (t_avg[j,i] - c_avg[j,i]))
    }
    mean_diff[i,1] <-  mean_diff[i,1] / total
  }
  # get pooled sd by using data after matching
  # sd_t <-numeric()
  # sd_c <-numeric()
  # sd_p <-numeric()
  # for(i in 1:ncol(X_mat)) {
  #   sd_t <- c(sd_t, sd(t_avg[,i]))
  #   sd_c <- c(sd_c, sd(c_avg[,i]))
  #   sd_p <- c(sd_p, sqrt((sd(t_avg[,i])^2 + sd(c_avg[,i])^2))/2)
  # }
  
  # get pooled sd by using data before matching
  X_mat_t = X_mat[t_ind, ]
  X_mat_c_before = X_mat[!t_ind, ]
  X_mat_t_var = apply(X_mat_t, 2, var)
  X_mat_c_before_var = apply(X_mat_c_before, 2, var)
  sd_p <- sqrt((X_mat_t_var + X_mat_c_before_var)/2)
  mean_diff[,2] <-  mean_diff[,1] / sd_p
  rownames(mean_diff) <- colnames(X_mat)
  colnames(mean_diff) <- c('mean difference', 'smd')
  return(list(obj_total = result$obj,  matched = strata_mip, num_of_stratum = num_of_stratum, smd_after_matching = mean_diff, time = time,
              constraint = Amat))
}
# mean_diff must be outcome of full_w_constraint
smd_plot <- function(X_mat, t_id, mean_diff, v_line = 0.1, legend_position = "topright") {
  X_mat_t = X_mat[t_id, ]
  X_mat_c_before = X_mat[-t_id, ]
  X_mat_c_before_mean = apply(X_mat_c_before, 2, mean)
  X_mat_t_mean = apply(X_mat_t, 2, mean)
  X_mat_t_var = apply(X_mat_t, 2, var)
  X_mat_c_before_var = apply(X_mat_c_before, 2, var)
  std_dif_before = (X_mat_t_mean - X_mat_c_before_mean)/sqrt((X_mat_t_var + X_mat_c_before_var)/2)
  abs_std_dif_before = abs(std_dif_before)
  n_aux = length(abs_std_dif_before)
  abs_std_dif_after = abs(mean_diff[,2])
  
  dotchart(abs_std_dif_before[n_aux:1], labels = colnames(X_mat)[n_aux:1], 
           cex = 0.7, pch = "", col = , main = "", xlim = c(0, 1), 
           xlab = "Absolute standardized differences in means")
  points(abs_std_dif_before[n_aux:1], y = 1:ncol(X_mat), cex = 0.9, 
         pch = 0)
  points(abs_std_dif_after[n_aux:1], y = 1:ncol(X_mat), cex = 0.8, 
         pch = 8, col = "blue")
  legend(legend_position, c("Before matching", "After matching"), 
         cex = 0.8, bty = "n", pch = c(0, 8), col = c("black", 
                                                      "blue"))
  abline(v = v_line, lty = 2)
} 

# example: COPD data
# Basic parameter
data <- read_csv('./copd.csv', ) # data with cut off = 0.15 
dist_mat = distmat(t_ind, X_mat) # 26 X 52
t_ind = data$treated
t_id = which(data$treated ==1)
data <- arrange(data, desc(treated)) # the data should be sorted in decreasing order by the treatment indicator
X_mat = as.matrix(data[ ,c('age', 'bmi', 'smk', 'py', 'Emph')])

mom_covs = X_mat
# mom_tols = cbind(2, 0.3, 0.05, 2, 0.001) # tolerance from 1-1 match
mom_tols = cbind(0.1, 0.1, 0.05, 5, 0.001)
w_constraints <- TRUE

matched <- full_w_constraint(X_mat = X_mat, dist_mat = dist_mat, t_ind = t_ind, max_control = 8,
                             mom_covs = X_mat, mom_tols = mom_tols, w_constraints = w_constraints)
matched$smd_after_matching
# total distance
matched$obj_total
# matched strata
matched$matched
# smd
matched$smd_after_matching
# time consume for solving MIP
matched$time
# plot covariate balance
smd_plot(X_mat = X_mat, t_id = t_id, mean_diff = matched$smd_after_matching, v_line = 0.02)
# analysis
data <- as.data.frame(data)
rownames(data) <- c(paste0("t", 1:sum(t_ind)), 
                    paste0("c", 1:sum(!t_ind)))
total_num <- nrow(data)
te <- 0
for(i in 1:matched$num_of_stratum) {
  size <- length(matched$matched[[i]])
  df <- data[matched$matched[[i]],]
  t_avg <- mean(df[startsWith(rownames(df), "t"), ]$Out_fSAD)
  c_avg <- mean(df[startsWith(rownames(df), "c"), ]$Out_fSAD)
  te <- te + ((t_avg - c_avg) * size)/total_num
}
te

######### compare with optimal full matching
full_out = matchit(treated ~ age + bmi + smk + py + Emph,
                   data = data, method = "full",
                   distance = dist_mat)
summary(full_out)
# compare minimized total distance, matched sample and CB
benchmark <- match.data(full_out)[c('subclass', 'treated')]
fullgrouped <- benchmark %>% mutate(id = row_number()) %>% group_by(subclass) %>% summarise(id = list(id))
strata_opt <- fullgrouped$id
total_dist_opt <- lapply(strata_opt, function(sublist) {
  i <- sublist[sublist <= 26]
  j <- sublist[sublist > 26]
  sum = 0
  for(x in i){
    for(y in j){
      sum <- sum + dist_mat[x, y-26]
    }
  }
  return(sum)
})

total_dist_opt  = sum(unlist(total_dist_opt))
total_dist_opt

# solved by mip wo const
matched_wo_const <- full_w_constraint(X_mat = X_mat, dist_mat = dist_mat, t_ind = t_ind, 
                                      mom_covs = X_mat, mom_tols = mom_tols, w_constraints = FALSE)

matched_wo_const$obj_total # same w/ min distance of optimal full matching
smd_plot(X_mat = X_mat, t_id = t_id, mean_diff = matched_wo_const$smd_after_matching, v_line = 0.02)
