library(igraph)
library(dplyr)

plot_network <- function(df, alpha = 0.05) {
  edges <- df %>%
    filter(P_Value < alpha) %>%
    select(Predictor, Effect)
  
  graph <- graph_from_data_frame(edges, directed = TRUE)
  
  plot(graph, 
       #edge.label = format(df$P_Value[df$P_Value < alpha], scientific = TRUE),
       edge.label.cex = 0.8,
       main = "Causal Diagram SDG",
       vertex.label.dist = 2,
       vertex.size = 10,
       vertex.label.cex = 0.8,
       vertex.color = "lightblue",
       edge.color = "gray",
       edge.arrow.size = 0.5)
}
plot_network(res)