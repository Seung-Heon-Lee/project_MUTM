rm(list = ls())
library('mipmatch')
data.full <- read.csv('/Users/lsh/Workspace/CausalInference/tutorial/data_demo_PRM_210310.csv')
data.male = data.full[data.full$male==1, ]
table(data.male$smoking_status_t1)
data.sub = data.male[data.male$smoking_status_t1 < 3, ]
#data.sub = data.sub[data.sub$GOLD_t1 < 3, ]
data.sub = data.sub[data.sub$Initial_Whole_Whole_Emph < 0.1, ]

data.sub$time15less = (data.sub$time_diff_t2_t1 < 1.5)
data.sub$time1530 = (data.sub$time_diff_t2_t1 >= 1.5 & data.sub$time_diff_t2_t1 < 3)
data.sub$time30more = (data.sub$time_diff_t2_t1 >= 3)
library(dplyr)
data.match.sub = dplyr::select(data.sub, age_t1, bmi_t1, smoking_status_t1, smoking_py_t1, GOLD_t1,
                               Initial_Whole_Whole_fSAD, Initial_Whole_Whole_Emph, Followup_Whole_Whole_fSAD, Followup_Whole_Whole_Emph)

colnames(data.match.sub) = c("age", "bmi", "smk", "py", "GOLD", "fSAD", "Emph", "Out_fSAD", "Out_Emph")
data.match.sub$treated  = (data.match.sub$fSAD < 0.150)
sum(data.match.sub$treated==1)
sum(data.match.sub$treated==0)
data.match.sub$smk[data.match.sub$smk==2] = 0

data.match.sub
sorted = arrange(data.match.sub, desc(treated))
write_csv(sorted,'./copd.csv')
t_ind = sorted$treated
t_id = which(sorted$treated ==1)
X_mat = cbind(sorted$age, sorted$bmi, sorted$smk, sorted$py, sorted$Emph)
dist_mat = distmat(t_ind, X_mat) # 26 X 52

full_out = matchit(treated ~ age + bmi + smk + py + Emph,
                   data = sorted, method = "full",
                   distance = 'mahalanobis')

full_out
summary(full_out)

# full matching using MIP
# Basic parameter
number_of_treated_nodes <- sum(t_ind)
number_of_control_nodes <- sum(!t_ind)
total <- number_of_control_nodes + number_of_treated_nodes
max_control <- number_of_control_nodes
max_treat <- number_of_treated_nodes
matched_control <- number_of_control_nodes # it can be adjusted
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
# must consider order of V(g)
sply_source <- number_of_treated_nodes*max_control
sply_sink <- -matched_control
sply_over <- -(sply_source - matched_control)
V(g)$supply <- c(sply_source, rep(0, total), sply_over, sply_sink) # supply of each control and treated node is 0
names(V(g)$supply) <- V(g)$name

# let's solve MIP
library(Rcplex)
# Coefficeints of the objective function,  variables are flow of each edge
cvec <- cost
# Construct the constraint matrix, for each node outgoing edge is 1, incoming edge is -1. If there's no edge than 0
Amat <- matrix(0, nrow=length(V(g)), ncol=length(E(g)))
for(i in 1:length(V(g))) {
  node <- V(g)$name[i]
  outgoing_edges <- incident(g, node, mode="out")
  incoming_edges <- incident(g, node, mode="in")
  
  # For outgoing edges of ith node, subtract flow from node supply
  Amat[i, outgoing_edges] <- 1 
  
  # For incoming edges, add flow to node supply
  Amat[i, incoming_edges] <- -1
}

# Vector of RHS of constraints (i.e, supply for each node, 0 if not specified)
bvec <- V(g)$supply
# Bounds of variables
lb <- rep(0, length(E(g)))
ub <- c(as.vector(capacity))
# the optimization direction
objsense = c('min')
# The direction of the inequality in each constraint. 'E' stands for equal
sense = 'E'
# Types of variables, Integer
vtype <- "I"
# Solve the model
result <- Rcplex(cvec = cvec, Amat = Amat, bvec = bvec, Qmat = NULL,
                 lb = lb, ub = ub, control = list(),
                 objsense = objsense, sense = sense, vtype = vtype, n = 1)


# get solution of full matching problem
library(tidyverse)
print(result)
edge_ids <- as_ids(E(g))
optimized_flows <- result$xopt
df <- data.frame(edges = edge_ids, optimized_flows = optimized_flows)
df_filtered <- df %>%
  filter(grepl("^t\\d+\\|c\\d+$", edges) & df$optimized_flows == 1)

df_filtered <- df_filtered %>%
  separate(edges, into = c("treated", "control"), sep = "\\|")

# Group by treated node and combine connected control nodes into a list
df_grouped <- df_filtered %>%
  group_by(treated) %>%
  summarise(control = list(control)) %>% 
  group_by(control) %>% 
  summarise(treated = list(treated)) %>%
  mutate(stratum = map2(control, treated, ~ paste(c(.x, .y), collapse = " ")))
df_grouped
strata_mip <- lapply(df_grouped$stratum, function(x) unlist(strsplit(x, split = " ")))
strata_mip
total_unit <- sum(unlist(lapply(strata_mip, length))) # = 78


# compare minimized total distance, matched sample and CB
benchmark <- match.data(full_out)[c('subclass', 'treated')]
fullgrouped <- benchmark %>% mutate(id = row_number()) %>% group_by(subclass) %>% summarise(id = list(id))
strata_opt <- fullgrouped$id
dist_opt <- lapply(strata_opt, function(sublist) {
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

total_dist_opt  = sum(unlist(total_dist_opt)) # 15.7
total_dist_opt

extract_number <- function(x) as.integer(sub("^[tc]", "", x))
pairs <- lapply(strata_mip, function(sublist) {
  t_values <- sapply(sublist[startsWith(sublist, "t")], extract_number)
  c_values <- sapply(sublist[startsWith(sublist, "c")], extract_number)
  
  expand.grid(t_values, c_values)
})

dist_mip <- sapply(pairs, function(x) {
  sapply(seq_len(nrow(x)), function(i) dist_mat[x[i, 1], x[i, 2]])
})

total_dist_mip = sum(unlist(dist_mip))
# minimized total distance in optmatch is 15.7
# minimized total distance in MIP is 12.6
matchit
