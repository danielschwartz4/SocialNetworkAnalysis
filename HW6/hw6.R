library(xUCINET)
library(sna)
library(igraph)
library(purrr)
network <- Corporation$NetworkInfo
attributes <- Corporation$Attributes
corperation <- Corporation$Corporation

undirected <- xSymmetrize(corperation, Type = "Max")

# 1.
cliques <- xCliquesMembership(undirected, Min = 6) 
cliques

# 2.
coCliques <- xCliquesCoMembership(undirected, Min = 6)

cluster <-
  xHierarchicalClustering(
    NET1 = coCliques,
    Input = "Similarities",
    Method = "average",
    TitleDendrogram = "Clique Overlap",
    NOC = 1:nrow(coCliques)
  )

graph <-
  graph_from_adjacency_matrix(adjmatrix = undirected,
                              mode = "undirected")

cluster_df <- data.frame(cluster)

modularityScores <-
  map_dbl(cluster_df, ~ modularity(x = graph,
                                      membership = .x))

sort(modularityScores, decreasing = TRUE)

#  3-cluster clustering:
split(x = rownames(cluster), f = cluster[, 2])

xDensity(undirected,
         ROWS = cluster[, 2],
         COLS = cluster[, 2])

# 4. 

GirvanNewman <- xGirvanNewman(NET1 = undirected)
split(x = rownames(GirvanNewman), f = GirvanNewman[, 1])
#  Note that SM_GirvanNewman has only one column

modularity(x = graph, membership = GirvanNewman[, "CL_4"])

xDensity(undirected,
         ROWS = GirvanNewman[, "CL_4"],
         COLS = GirvanNewman[, "CL_4"])

table(GirvanNewman[, "CL_4"], attributes$office)
table(GirvanNewman[, "CL_4"], attributes$projects)
table(GirvanNewman[, "CL_4"], attributes$seniority)

tapply(X = attributes$office, INDEX = GirvanNewman[, "CL_4"], FUN = mean)
tapply(X = attributes$projects, INDEX = GirvanNewman[, "CL_4"], FUN = mean)
tapply(X = attributes$seniority, INDEX = GirvanNewman[, "CL_4"], FUN = mean)






