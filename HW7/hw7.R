library(xUCINET)
library(sna)

load("~/documents/college/socialnetworkanalysis/socialnetworkanalysisrepo/hw7/Movies.Rdata")

collab <- Movies$Collab
producer <- Movies$Attributes$Producer
composer <- Movies$Attributes$composer

Collab_Bipartite <-
  xTwoModeToBipartite(NET2 = collab,
                      First = "RowNodes")

# 1.
xDensity(Collab_Bipartite)
dim(collab)
xComponents(collab)

# 2. 

ProducerProj <-
  xTwoModeToOneMode(collab,
                    Measure = "CrossProd")

ComposerProj <-
  xTwoModeToOneMode(t(collab),
                    Measure = "CrossProd")

ComposerProj[1:10, 1:10]
# Max is 5
which(collab[, "C17"] * collab[, "C31"] == 1)


# 3. 
# -- Degree
ProducerDegree <-
  xDegreeCentrality(collab)
ComposerDegree <-
  xDegreeCentrality(t(collab))

ProducerDegree[1:10, ]
ComposerDegree[1:10, ]

head(sort(ProducerDegree[, 1], decreasing = TRUE))
head(sort(ComposerDegree[, 1], decreasing = TRUE))

# -- Betweeness

Collab_Betweenness <- xBetweennessCentrality(Collab_Bipartite)
Collab_Betweenness_adj <- Collab_Betweenness[, 1] / 2

Collab_Betweenness_Producers <-
  Collab_Betweenness_adj[1:NROW(collab)] # 1 to N elements
Collab_Betweenness_Composers <-
  Collab_Betweenness_adj[(NROW(collab) + 1):NROW(Collab_Bipartite)] # N + 1 to N + M elements

Collab_Betweenness_Producers[1:10]
Collab_Betweenness_Composers[1:10]

head(sort(Collab_Betweenness_Producers, decreasing = TRUE))
head(sort(Collab_Betweenness_Composers, decreasing = TRUE))


# 4. 
par(mar = c(0, 0, 0, 0))
set.seed(02138)
Collab_PlotIncidence <- gplot(
  collab,
  gmode = "twomode",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80"
)

Collab_MostCentralBetwenness_Producers <-
  names(sort(Collab_Betweenness_Producers, decreasing = TRUE)[1:2])

Collab_MostCentralBetwenness_Composers <-
  names(sort(Collab_Betweenness_Composers, decreasing = TRUE)[1:2])

Collab_color <- rep("red", NROW(Collab_Bipartite))

Collab_MostCentralBetwenness_ProducerIndex <-
  match(Collab_MostCentralBetwenness_Producers,
        rownames(Collab_Bipartite))

Collab_color[Collab_MostCentralBetwenness_ProducerIndex] <- "green"

#  Repeat the steps above for the most central events to color them blue:

Collab_MostCentralBetwenness_ComposerIndex <-
  match(Collab_MostCentralBetwenness_Composers,
        rownames(Collab_Bipartite))

Collab_color[Collab_MostCentralBetwenness_ComposerIndex] <- "blue"

set.seed(02138)
gplot(
  collab,
  gmode = "twomode",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80",
  vertex.col = Collab_color,
)

# 5.
Collab_CP <- xDualCorePeriphery(collab)
Collab_CP

#  Let's check the membership of women in the core/periphery:

split(rownames(collab), Collab_CP[[1]])

#  Similarly for composers:

split(colnames(collab), Collab_CP[[2]])

xDensity(collab,
         ROWS = Collab_CP[[1]],
         COLS = Collab_CP[[2]])

library(blockmodeling)
plot.mat(collab,
         clu = list(Collab_CP[[1]],
                    Collab_CP[[2]]))
Collab_color2 <- rep("red", NROW(Collab_Bipartite))
Collab_color2[match(names(which(Collab_CP[[1]] == 1)), rownames(Collab_Bipartite))] <- "green"
Collab_color2[match(names(which(Collab_CP[[2]] == 1)), rownames(Collab_Bipartite))] <- "blue"

set.seed(02138)
gplot(
  collab,
  gmode = "twomode",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80",
  vertex.col = Collab_color2,
)





