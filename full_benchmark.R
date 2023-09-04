library(optmatch)
library(MatchIt)
library(mipmatch)
library(igraph)
lalonde # 185 treated
lalonde$black = lalonde$race == "black"
lalonde$hisp =lalonde$race == "hispan"
data = lalonde
t_ind = data$treat
t_id = which(t_ind==1)
X_mat = cbind(data$age, data$educ, data$married, data$nodegree, data$re74, data$re75, data$black, data$hisp)
full_out = matchit(treat ~ age + educ + married + black + hisp + nodegree + re74 +re75,
                   data = lalonde, method = "full",
                   distance = 'mahalanobis', max_controls = 30, max_treat = 30)

# full_out = matchit(treat ~ age + educ + married + black + hisp + nodegree + re74 +re75,
#                    data = lalonde, method = "full",
#                    distance = dist_mat, max_controls = 30)
                   

summary(full_out)

# full matching using MIP
# Basic parameter
number_of_treated_nodes <- sum(t_ind)
number_of_control_nodes <- sum(!t_ind)
total <- number_of_control_nodes + number_of_treated_nodes
max_control <- 30
max_treat <- 1000
matched_control <- number_of_control_nodes # it can be adjusted
dist_mat = distmat(t_ind, X_mat) # 185 X 429
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
sply_sink <- -matched_control
sply_over <- -(sply_source - matched_control)
V(g)$supply <- c(sply_source, rep(0, total), sply_sink, sply_over) # supply of each control and treated node is 0
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
library(dplyr)
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
  summarise(control = list(control)) 

# have to make matched stratum


