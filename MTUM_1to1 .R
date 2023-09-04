# Multi Threshold Update Match 
library(igraph) # to make graph
library(mipmatch) # use distmat function to construct distance matrix
library(Rcplex) # to solve IP
library(dplyr) # for data wrangling
library(tidyverse)
library("lattice") # for plotting
distmat <- function (z, X, digits = 1) {
  X = as.matrix(X)
  n = dim(X)[1]
  rownames(X) = 1:n
  k = dim(X)[2]
  m = sum(z)
  for (j in 1:k) X[, j] = rank(X[, j])
  cv = cov(X)
  vuntied = var(1:n)
  diag(cv)[diag(cv) == 0] = 0.01
  rat = sqrt(vuntied/diag(cv))
  cv = diag(rat) %*% cv %*% diag(rat)
  out = matrix(NA, m, n - m)
  Xc = X[z == 0, ]
  Xt = X[z == 1, ]
  rownames(out) = rownames(X)[z == 1]
  colnames(out) = rownames(X)[z == 0]
  library(MASS)
  icov = ginv(cv)
  if(m == 1) {out = mahalanobis(Xc, Xt, icov, inverted = TRUE)}
  else {
    for (i in 1:m) {
      out[i, ] = mahalanobis(Xc, Xt[i, ], icov, inverted = TRUE)
    } 
  }
  dist_mat = out
  dist_mat = dist_mat/mean(dist_mat)
  dist_mat = round(dist_mat/mean(dist_mat), digits)
  dist_mat
}
# for COPD data
# Basic parameter
data <- as.data.frame(read_csv('/Users/lsh/Workspace/CausalInference/project/MIP_match/data/copd.csv')) # cutoff: 0.15, decendently sorted by treated
sorted_fSAD <- sort(unique(data$fSAD))
# starts with 3 treated unit
data$treated <- data$fSAD <= sorted_fSAD[3]
# data$treated <- data$fSAD <= sorted_fSAD[27]
data <- data %>% arrange(desc(treated)) %>% arrange(fSAD) %>% mutate(unit_id = seq_len(nrow(data)))
t_ind = data$treated
t_id = which(data$treated ==1)
num_t = sum(t_ind)
X_mat = as.matrix(data[ ,c('age', 'bmi', 'smk', 'py', 'Emph')])
dist_mat = distmat(t_ind, X_mat)
mom_covs = X_mat
mom_tols = cbind(1, 0.3, 0.05, 2, 0.001) 
X_mat_t = X_mat[t_ind, ]
X_mat_c_before = X_mat[!t_ind, ]
X_mat_t_var = apply(X_mat_t, 2, var)
X_mat_c_before_var = apply(X_mat_c_before, 2, var)
sd_p <- sqrt((X_mat_t_var + X_mat_c_before_var)/2)
tol_1 <- sd_p * 0.1
tol_2 <- sd_p * 0.2
tol_1
tol_2

# initial 1-1 matching with constraints, starts with 3 treated unit, to achieve CB
mom_tols = tol_2
matched <- mipmatch::allmatch(dist_mat = dist_mat, t_ind = t_ind, n_matches = 1 , mom_covs = mom_covs, mom_tols = mom_tols)
jotplot(X_mat, t_id, matched$c_id, v_line = 0.1, legend_position = "topright")
meantab(X_mat, t_ind, t_id, matched$c_id, digits = 2)
matched$c_id
matched_control <- data[matched$c_id,] %>%
  mutate(matched_treated = data$unit_id[1:3])  %>%
  mutate(dist = dist_mat[cbind(matched_treated, unit_id-num_t)])
unmatched_control <- data[-matched$c_id,] %>% filter(treated==FALSE)
matched_control
unmatched_control

# change one control unit to treated
result <- list()

#
result[[3]] <- matched_control
for (i in seq_len(length(sorted_fSAD))) {
  if(i>3 && i<=39) {
    moved <- data
    moved$treated <- data$fSAD <= sorted_fSAD[i]
    moved <- moved %>% arrange(desc(treated))
    moved_unit <- moved %>% filter(treated == TRUE) %>% arrange(desc(fSAD)) %>% head(1)
    t_ind <- moved$treated
    # case1: changed unit was matched control
    if(nrow(anti_join(matched_control, moved_unit, by = "unit_id")) != nrow(matched_control)) {
      print('case1')
      released_treated_id <- unlist(matched_control %>% filter(unit_id == moved_unit$unit_id) %>% select(matched_treated))
      matched_control <- anti_join(matched_control, moved_unit, by = "unit_id")
      partial_treated <- rbind(moved[moved$unit_id == released_treated_id,], moved_unit)
      partial <- rbind(partial_treated, unmatched_control)
      partial_X <- as.matrix(partial[ ,c('age', 'bmi', 'smk', 'py', 'Emph')])
      partial_t_ind <- partial$treated
      partial_dist <- distmat(partial_t_ind, partial_X)
      partial_matched <- allmatch(partial_dist, partial_t_ind, n_matches = 1)
      newly_matched <- partial[partial_matched$c_id,] %>%
        mutate(matched_treated = c(released_treated_id, moved_unit$unit_id)) %>% 
        mutate(dist = partial_dist[cbind(c(1,2), partial_matched$c_id - 2)])
      matched_control <- rbind(matched_control, newly_matched)
      unmatched_control <- moved %>% filter(!unit_id %in% matched_control$unit_id) %>%
        filter(treated==FALSE)
    }
  
    # case2: changed unit was unmatched control
    if(nrow(anti_join(unmatched_control, moved_unit, by = "unit_id")) != nrow(unmatched_control)) {
      print('case2')
      unmatched_control <- anti_join(unmatched_control, moved_unit, by = "unit_id")
      partial <- rbind(moved_unit, unmatched_control)
      partial_X <- as.matrix(partial[ ,c('age', 'bmi', 'smk', 'py', 'Emph')])
      partial_t_ind <- partial$treated
      partial_dist <- distmat(partial_t_ind, partial_X)
      partial_matched <- allmatch(partial_dist, partial_t_ind, n_matches = 1)
      newly_matched <- partial[partial_matched$c_id,] %>%
        mutate(matched_treated = moved_unit$unit_id) %>%
        mutate(dist = partial_dist[partial_matched$c_id - 1])
      matched_control <- rbind(matched_control, newly_matched)
      unmatched_control <- moved %>% filter(!unit_id %in% matched_control$unit_id) %>%
        filter(treated==FALSE)
    }
    # matched control for each cut off
    result[[i]] <- matched_control
  }
}
# total distance
total_dist <- c()
for(i in 3:39) {
  df <- result[[i]]
  total_dist <- c(total_dist, unlist(sum(df$dist)))
}


# smd
smd <- data.frame(X_mat[3:39,])
for(i in 3:39) {
  df <- result[[i]]
  smd[i-2,] <- colMeans(X_mat[df$matched_treated,] - X_mat[df$unit_id,]) / sd_p
}
smd # treated 3개에서 39개까지 CB 변화

# treatment effect
te <- c()
for(i in 3:39) {
  df <- result[[i]]
  treat <- data %>% filter(data$unit_id %in% df$matched_treated) %>% select(Out_Emph)
  control <- data %>% filter(data$unit_id %in% df$unit_id) %>% select(Out_Emph)
  te <- c(te, colMeans(treat - control))
}
te

# benchmark, solve mip for each cutoff
# 각 cutoff에서 어떻게 tolerance를 줘야 알고리즘과 알맞게 비교할 수 있을까?
benchmark <- list()
moved <- data
mom_tols <- cbind(1, 0.3, 0.05, 2, 0.001)
# three benchmark option
which = 3
for (i in seq_len(length(sorted_fSAD))) {
  moved$treated  = (data$fSAD <= sorted_fSAD[i])
  if(sum(moved$treated==1) >= 3 && sum(moved$treated==1) <= sum(moved$treated==0) ) { # if only one treated or control, no need to match
    sorted = arrange(moved, desc(treated))
    t_ind = sorted$treated
    t_id = sorted$unit_id[sorted$treated==1]
    X_mat = as.matrix(sorted[ ,c('age', 'bmi', 'smk', 'py', 'Emph')])
    dist_mat = distmat(t_ind, X_mat)
    if(which==1) {result = allmatch(dist_mat, t_ind, n_matches = 1)}
    if(which==2){
      if(i <=26) {
        result = allmatch(dist_mat, t_ind, n_matches = 1,
                          mom_covs = mom_covs, mom_tols = sd_p * 0.11,
                          enforce_constraints = TRUE)
      }
      else {result = allmatch(dist_mat, t_ind, n_matches = 1)}
    }
    if(which == 3) {
      if(i <=26) {
        result = allmatch(dist_mat, t_ind, n_matches = 1,
                          mom_covs = mom_covs, mom_tols = sd_p * 0.11,
                          enforce_constraints = TRUE)
      }
      else if (i <= 28) {
        result = allmatch(dist_mat, t_ind, n_matches = 1,
                          mom_covs = mom_covs, mom_tols = sd_p * 0.2,
                          enforce_constraints = TRUE)
      }
      else if (i <= 31) {
        result = allmatch(dist_mat, t_ind, n_matches = 1,
                          mom_covs = mom_covs, mom_tols = sd_p * 0.3,
                          enforce_constraints = TRUE)
      }
      else {result = allmatch(dist_mat, t_ind, n_matches = 1)} # treated가 32개 넘으면 w/o constraint
    }
    c_id = sorted[result$c_id, 'unit_id'] #sorted의 row기준 index를 unit_id로 바꿔줘야함
    benchmark[[i]] <- list(t_id = t_id, c_id = c_id, dist = result$obj_total)
  }
}

# t_id, c_id 모두 unit_id 기준, dist_mat는 sorted index 기준, distance를 알고싶으면 sorted의 row index로 바꿔서
# matching을 매번 새로 찾았을 때 기준
benchmark_df <- as.data.frame(do.call(rbind, benchmark))
benchmark_df$t_id
benchmark_df$c_id
benchmark_df$dist
total_dist_benchmark <- unlist(benchmark_df$dist)
smd_benchmark <- data.frame(X_mat[3:39,])
for(i in 3:39) {
  df <- benchmark[[i]]
  smd_benchmark[i-2,] <- colMeans(X_mat[df$t_id,] - X_mat[df$c_id,]) / sd_p
}
te_benchmark <- c()
for(i in 3:39) {
  df <- benchmark[[i]]
  treat <- data %>% filter(data$unit_id %in% df$t_id) %>% select(Out_Emph)
  control <- data %>% filter(data$unit_id %in% df$c_id) %>% select(Out_Emph)
  te_benchmark <- c(te_benchmark, colMeans(treat-control))
}
te_benchmark

# visualize
library(ggplot2)
library(gridExtra)
t_count <- 3:39
# total distance
plot(t_count, total_dist, type="l", col="blue", lwd=2, ylab="total_dist", xlab="num of t", main="compare total distance")
lines(t_count, total_dist_benchmark, type="l", col="red", lwd=2)
legend("topright", legend=c("algorithm", "benchmark"), col=c("blue", "red"), lwd=2)

# smd
plots <- lapply(names(smd), function(colname) {
  ggplot() +
    geom_line(data=smd, aes_string(x=t_count, y=paste("abs(", colname, ")")), color="blue") +
    geom_line(data=smd_benchmark, aes_string(x=t_count, y=paste("abs(", colname, ")")), color="red") +
    labs(title=paste("Comparison of smd", colname), x="num of t", y="Absolute Value") +
    geom_hline(aes(yintercept=0.1), color="black", linetype="dashed") +
    theme_minimal()
})
do.call(grid.arrange, plots)

# treatment effect
plot(t_count, te, type="l", col="blue", lwd=1, ylab="treatment effect", xlab="num of t", main="compare treatment effect",
     ylim = c(-0.03,0))
lines(t_count, te_benchmark, type="l", col="red", lwd=1)
legend("topright", legend=c("algorithm", "benchmark"), col=c("blue", "red"), lwd=2)

sorted_fSAD[27]
