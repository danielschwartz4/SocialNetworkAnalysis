library(xUCINET)
library(sna)


attributes <- Zachary_KarateClub$Attributes
network <- Zachary_KarateClub$Connection


# Degree centrality
dcentrality <- xDegreeCentrality(network)

# Order
dcentralitydf <- data.frame(dcentrality)
dcentralitydf[order(dcentralitydf$Degree),]

dmean <- mean(dcentrality)

# dc vis
par(mar= c(1, 1, 1, 1))
set.seed(02138)
gplot(network, 
      displaylabels = TRUE,
      vertex.cex = dcentrality / dmean,
      vertex.col = attributes$Club + 1)

# Betweeness centrality
bcentrality <- xBetweennessCentrality(network)
bcentralitydf
bcentralitydf <- data.frame(bcentrality)
bcentralitydf[order(bcentralitydf$Betweenness),]

# Eigenvector centrality
ecentrality <- xEigenvectorCentrality(network)
ecentrality
emean <- mean(ecentrality)

# ec vis
# Vis
gplot(network, 
      displaylabels = TRUE,
      vertex.cex = ecentrality * 10,
      vertex.col = attributes$Club + 1)

# Closeness centrality (default freeman, normalize true)
# ccentrality <- centralization(network,  
#                               closeness,
#                               mode="graph")
ccentrality <- xClosenessCentrality(network)


# Aggregate in one table
cent_scores <- data.frame(dcentrality[,1], 
                          ccentrality[,1],
                          bcentrality[,1],
                          ecentrality[,1])
colnames(cent_scores)<-c("Degree","Closeness", "Betweenness","Eigenvector")
round(cor(cent_scores, use = "complete.obs"), digits = 3)


# degree
centralization(network, degree, normalize = TRUE)
centralization(network, degree, normalize = FALSE)
# betweenness
centralization(network, betweenness, normalize = TRUE)
centralization(network, betweenness, normalize = FALSE)








