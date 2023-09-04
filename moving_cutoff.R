detach()
rm(list = ls())
library('mipmatch')
data.full <- read.csv('/Users/lsh/Workspace/CausalInference/tutorial/data_demo_PRM_210310.csv')
data.male = data.full[data.full$male==1, ]
table(data.male$smoking_status_t1)
data.sub = data.male[ddata.male$smoking_status_t1 < 3, ]
#data.sub = data.sub[data.sub$GOLD_t1 < 3, ]
data.sub = data.sub[data.sub$Initial_Whole_Whole_Emph < 0.1, ]

data.sub$time15less = (data.sub$time_diff_t2_t1 < 1.5)
data.sub$time1530 = (data.sub$time_diff_t2_t1 >= 1.5 & data.sub$time_diff_t2_t1 < 3)
data.sub$time30more = (data.sub$time_diff_t2_t1 >= 3)
library(dplyr)
data.match.sub = dplyr::select(data.sub, age_t1, bmi_t1, smoking_status_t1, smoking_py_t1, GOLD_t1,
                               Initial_Whole_Whole_fSAD, Initial_Whole_Whole_Emph, Followup_Whole_Whole_fSAD, Followup_Whole_Whole_Emph)

colnames(data.match.sub) = c("age", "bmi", "smk", "py", "GOLD", "fSAD", "Emph", "Out_fSAD", "Out_Emph")
data.match.sub$smk[data.match.sub$smk==2] = 0
data.match.sub$unit_id <- seq_len(nrow(data.match.sub))
cutoff <- 0.150
data.match.sub$treated  = (data.match.sub$fSAD < cutoff)
sum(data.match.sub$treated==1)
sum(data.match.sub$treated==0)

sorted = arrange(data.match.sub, desc(treated))
t_ind = sorted$treated
t_id = which(sorted$treated ==1)
X_mat = cbind(sorted$age, sorted$bmi, sorted$smk, sorted$py, sorted$Emph)
dist_mat = distmat(t_ind, X_mat)

out_wo_constraints = allmatch(dist_mat, t_ind, n_matches = 1)
meantab(X_mat,t_ind,t_id,out_wo_constraints$c_id,digits=2)

# use tolerance, cutoff = 0.15
mom_covs = cbind(age, bmi, smk, py, Emph)
mom_tols = cbind(1, 0.3, 0.05, 2, 0.001)
out_w_constraints = allmatch(dist_mat, t_ind, n_matches = 1,
                             mom_covs = mom_covs, mom_tols = mom_tols,
                             enforce_constraints = TRUE)
meantab(X_mat, t_ind, t_id, out_w_constraints$c_id, digits = 2)

c_matched_w_const = out_w_constraints$c_id

mipmatch::jotplot(X_mat, t_id, c_matched_w_const, v_line = 0.1)
mipmatch::jotplot(X_mat, t_id, out_wo_constraints$c_id, v_line = 0.1)
mean(Emph)
sd(Emph)

# use tolerance, cutoff = 0.15
cutoff <- 0.155
data.match.sub$treated  = (data.match.sub$fSAD < cutoff)
sum(data.match.sub$treated==1)
sum(data.match.sub$treated==0)

sorted = arrange(data.match.sub, desc(treated))
t_ind = sorted$treated
t_id = which(sorted$treated ==1)
X_mat = cbind(sorted$age, sorted$bmi, sorted$smk, sorted$py, sorted$Emph)
dist_mat = distmat(t_ind, X_mat)

# as cutoff increase, a control unit changed to treated unit. find each point(cutoff) that adjust composition of treated and control set
fSAD <- as.vector(data.match.sub$fSAD)
sorted_fSAD <- sort(unique(fSAD))

change_points <- c()
prev_treated <- NULL

for (cut_off in sorted_fSAD) {
  treated <- sum(fSAD <= cut_off)
  if (is.null(prev_treated) || prev_treated != treated) {
    change_points <- c(change_points, cut_off)
    prev_treated <- treated
  }
}
change_points
data.match.sub
mom_covs = cbind('age', 'bmi', 'smk', 'py', 'Emph')
mom_tols = cbind(1, 0.3, 0.05, 2, 0.001) # 0.15 기준으로 balance를 잘 맞출 수 있었던 tolerance
matched <- list()
for (cut_off in change_points) {
  data.match.sub$treated  = (data.match.sub$fSAD < cut_off)
  if(sum(data.match.sub$treated==1) > 1 && sum(data.match.sub$treated==1) <= sum(data.match.sub$treated==0) ) { # if only one treated or control, no need to match
    sorted = arrange(data.match.sub, desc(treated))
    t_ind = sorted$treated
    t_id = sorted$unit_id[sorted$treated==1]
    X_mat = cbind(sorted$age, sorted$bmi, sorted$smk, sorted$py, sorted$Emph)
    dist_mat = distmat(t_ind, X_mat)
    result = allmatch(dist_mat, t_ind, n_matches = 1,
                      mom_covs = X_mat, mom_tols = mom_tols,
                      enforce_constraints = TRUE)
    c_id = sorted[result$c_id, 'unit_id'] #sorted의 row기준 index를 unit_id로 바꿔줘야함
    matched[[as.character(cut_off)]] <- list(t_id = t_id, c_id = c_id, dist = result$obj_total)
  }
}
# t_id, c_id 모두 unit_id 기준, dist_mat는 sorted index 기준, distance를 알고싶으면 sorted의 row index로 바꿔서
# matching을 매번 새로 찾았을 때 기준
matched_df <- as.data.frame(do.call(rbind, matched))
matched_df$t_id
matched_df$c_id
matched_df$dist

