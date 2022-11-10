#  Sociology 2275: Social Network Analysis
#  Fall 2022
#  Section 7 - Date: November 4, 2022
#
#  Course Instructor:  Peter V. Marsden <pvm@wjh.harvard.edu>
#
#  Teaching Fellow:    Derick S. Baum <derick_baum@g.harvard.edu>,
#                      Office hours: Thursdays 3 PM - 5 PM, or by appointment
#





#  1. Setting your working directory and preparing data ---------------------------------------------

#  Follow the usual instructions
#  to set your working directory to the folder where the R script is saved.
#  We will use the rstudioapi package to accomplish this:

install.packages("rstudioapi") # Install the package if you haven't done so already

#  Load the package:

library(rstudioapi)

#  Use Cmd + Enter to run the line below. This will get the path to the current R script:

current_path <- getActiveDocumentContext()$path

#  Now set the working directory to the folder where the script is saved.
#  setwd() sets the working directory.
#  dirname() gets the directory name.

setwd(dirname(current_path))

#  Let's check if we were successful:

getwd() # This should be the location of your current R script

#  Today we'll study two-mode data using a famous example:
#  The Davis et al. (1941) data on the attendance of 18 Southern women at 14 society events,
#  discussed at length in Chapter 13 of the Borgatti et al. textbook

library(xUCINET)


#  2. Basic descriptive statistics ---------------------------------------------

## Incidence matrix and bipartite network ---------------------------------------------

#  First store the incidence matrix and the bipartite network in separate objects:

SouthernWomen_Incidence <- Davis_SouthernWomen$Attendance

#  The xTwoModeToBipartite() function converts the incidence matrix into a one-mode network
#  that includes one row/column for each row and column of the two-mode network.
#  So if there are N rows and M columns in the incidence matrix, 
#  this "bipartite" matrix will have N+M rows and columns. 
#  In this network, there are no within-type edges: 
#  all relationships among the N row nodes are defined as 0, 
#  and likewise all relationships among the M column nodes are defined as 0:

SouthernWomen_Bipartite <-
  xTwoModeToBipartite(NET2 = SouthernWomen_Incidence,
                      First = "RowNodes")

#  NET2 is the incidence matrix
#  First indicates whether the "RowNodes" or the "ColumnNodes" should appear first in the bipartite network. 
#  It defaults to "RowNodes", meaning the row entities in the incidence matrix are in rows 1 to N of the bipartite network, 
#  and the column entities in the incidence matrix in rows N+1 to N+M of the bipartite network. 
#  With First set to "ColumnNodes", the column entities in the incidence matrix appear in rows 1 to M of the bipartite network, 
#  while the row entities in the incidence matrix appear in rows M+1 to N+M of the bipartite matrix


## Size ---------------------------------------------

#  The network size in the case of two-mode data is the number of row and column entities
#  This can be computed using the NROW() and NCOL() functions:

NROW(SouthernWomen_Incidence)
NCOL(SouthernWomen_Incidence)

#  Or directly using the dim() function (the number of dimensions of the matrix):

dim(SouthernWomen_Incidence)


## Density ---------------------------------------------

#  We can compute density directly by applying the xDensity() function to the incidence matrix:

xDensity(SouthernWomen_Incidence)

#  There are 89 edges in the network, as given by the sum of the entries in the incidence matrix: 

SouthernWomen_edges <- sum(SouthernWomen_Incidence)
SouthernWomen_edges

#  We can compute density by dividing this value by the the maximum possible,
#  the product of the number of row and column entities:

SouthernWomen_edges / (NROW(SouthernWomen_Incidence) * NCOL(SouthernWomen_Incidence))

#  Note that the value will be too small if we use the bipartite network instead:

xDensity(SouthernWomen_Bipartite)

#  This is because the function does not take the structural 0 values for within-type relationships into account


## Number of components ---------------------------------------------

#  The number of components can be calculated by applying the xComponents() function to the incidence matrix.
#  It gives the number of components, minimum/maximum component size, etc. for the projection matrix for row entities ("Mode A"), 
#  the projection matrix for column entities ("Mode B"), and the bipartite network ("Both")

xComponents(SouthernWomen_Incidence)

#  We see one component for Mode A, Mode B, and the bipartite network,
#  meaning every pair of women has attended at least one event in common,
#  and every pair of events shared the attendance of at least one woman

#  Let's check how the output above would have changed if DOROTHY had only attended event E01,
#  and this event had only been attended by DOROTHY.
#  For this I will create an alternative incidence matrix that imagines this scenario:

SouthernWomen_IncidenceAlt <- SouthernWomen_Incidence
SouthernWomen_IncidenceAlt[, 1] <- 0
SouthernWomen_IncidenceAlt["DOROTHY", ] <- 0
SouthernWomen_IncidenceAlt["DOROTHY", 1] <- 1
SouthernWomen_IncidenceAlt

#  Let's pass the alternative incidence matrix to xComponents():

xComponents(SouthernWomen_IncidenceAlt)

#  We can now see two components for Mode A (women) and Mode B (events).
#  DOROTHY is now isolated because she didn't appear together in any event with any other woman.
#  Similarly, event E01 is an isolate because it didn't share the attendance of anyone with any other event 

#  Note that the component distribution for the bipartite network indicates two components,
#  one with size 2 and one with size 30. The component of size 2 consists of DOROTHY and event E01,
#  while all other women and events are in the largest component

#  Let's visualize dichotomized versions of the one-mode projections for this alternative incidence matrix
#  (with entries equal to 1 if two women appear together in at least one event; 
#  and similarly, if a pair of events shares the attendance of at least one woman)

#  Let's first visualize the woman-by-woman projection:

library(sna)
SouthernWomen_WomenProjAlt <- xTwoModeToOneMode(SouthernWomen_IncidenceAlt)
par(mar = c(0, 0, 0, 0))
set.seed(02138)
gplot(
  SouthernWomen_WomenProjAlt >= 1,
  gmode = "graph",
  usearrows = FALSE,
  displaylabels = TRUE,
  edge.col = "gray80",
  label.cex = 0.8
)

#  And now the event-by-event projection:
SouthernWomen_EventProjAlt <- xTwoModeToOneMode(t(SouthernWomen_IncidenceAlt))

set.seed(02138)
gplot(
  SouthernWomen_EventProjAlt >= 1,
  gmode = "graph",
  usearrows = FALSE,
  displaylabels = TRUE,
  edge.col = "gray80",
  label.cex = 0.8
)

#  Finally, we visualize the bipartite network (more on visualization later; 
#  for now, just focus on the output):

set.seed(02138)
gplot(
  SouthernWomen_IncidenceAlt,
  gmode = "twomode",
  usearrows = FALSE,
  displaylabels = TRUE,
  edge.col = "gray80",
  label.cex = 0.8
)




#  2. One-mode projections ---------------------------------------------

## Constructing the one-mode projections from the two-mode data ---------------------------------------------

#  We use the xTwoModeToOneMode() function to construct the one-mode projections
#  First, we construct the woman-by-woman projection whose cells indicate
#  the number of times two women attended an event together

SouthernWomen_WomenProj <-
  xTwoModeToOneMode(SouthernWomen_Incidence,
                    Measure = "CrossProd")

#  Using the default value for the Measure argument,
#  this is equivalent to multiplying the incidence matrix by its transpose:

SouthernWomen_WomenProj_manual <-
  SouthernWomen_Incidence %*% t(SouthernWomen_Incidence)

#  %*% indicates matrix multiplication
#  Check if the results are identical:

identical(SouthernWomen_WomenProj,
          SouthernWomen_WomenProj_manual)

#  To create a one-mode projection linking the column entities with one another, 
#  we need to apply xTwoModeToOneMode() to the transpose of the incidence matrix:

SouthernWomen_EventProj <-
  xTwoModeToOneMode(t(SouthernWomen_Incidence),
                    Measure = "CrossProd")

#  This is equivalent to multiplying the transpose of the incidence matrix by the original matrix:

SouthernWomen_EventProj_manual <-
  t(SouthernWomen_Incidence) %*% SouthernWomen_Incidence

identical(SouthernWomen_EventProj, 
          SouthernWomen_EventProj_manual)

#  The entries in this matrix show how many women attended both events in common


## Questions about the event-by-event projection ---------------------------------------------

#  Consider the one-mode projection showing how many women attended a pair of events in common
#  Display its first 5 rows and columns:

SouthernWomen_EventProj[1:5, 1:5]

#  We can see that events E01 and E05 share the attendance of 3 women.
#  Let's find out who they are:

which(SouthernWomen_Incidence[, "E01"] * SouthernWomen_Incidence[, "E05"] == 1)

#  EVELYN, LAURA and BRENDA attended events E01 and E05 in common

#  Which pair(s) of events shared the attendance of most women?
#  Let's first transform the diagonal elements of the one-mode projection into NAs to facilitate computation
#  Note that the diagonal elements refer to how many women attended that event.
#  They are, therefore, equal to the column sums of the incidence matrix:

identical(diag(SouthernWomen_EventProj),
          colSums(SouthernWomen_Incidence))

#  To preserve the original one-mode projection, I will transform a copy of it:

SouthernWomen_EventProj_copy <- SouthernWomen_EventProj
diag(SouthernWomen_EventProj_copy) <- NA

#  To locate the row and column with the maximum value in this network,
#  we will use the which() function setting its arr.ind argument to TRUE.
#  This allows us to get the row and column indices of the maximum value.
#  Otherwise, the matrix is treated as one big vector and indices returned are not very helpful:

which(
  SouthernWomen_EventProj_copy == max(SouthernWomen_EventProj_copy, na.rm = TRUE),
  arr.ind = TRUE
)

#  We see that events E08 and E09 shared the attendance of most women (9):

max(SouthernWomen_EventProj_copy, na.rm = TRUE)
SouthernWomen_EventProj_copy[8, 9]

#  These women are:

which(SouthernWomen_Incidence[, "E08"] * SouthernWomen_Incidence[, "E09"] == 1)




#  3. Centrality scores ---------------------------------------------

## Degree centrality ---------------------------------------------

#  We can compute the raw and normalized degree centrality scores directly from the incidence matrix
#  By default, the function will give centrality scores for the row nodes:

SouthernWomen_DegreeWomen <-
  xDegreeCentrality(SouthernWomen_Incidence)

#  The raw scores equal the row sums of the incidence matrix:

identical(SouthernWomen_DegreeWomen[, 1],
          rowSums(SouthernWomen_Incidence))

#  The normalized scores equal the row sums divided by maximum number of events a woman could have attended
#  (this is simply the number of columns in the incidence matrix):

identical(
  SouthernWomen_DegreeWomen[, 2],
  rowSums(SouthernWomen_Incidence) / NCOL(SouthernWomen_Incidence)
)

#  To obtain the degree centrality scores for the events, 
#  we need to pass the transposed incidence matrix to xDegreeCentrality():

SouthernWomen_DegreeEvents <-
  xDegreeCentrality(t(SouthernWomen_Incidence))

#  Note that using the bipartite matrix to compute the degree centrality scores
#  yields the correct raw scores. However, the normalized ones will be too small
#  because the function does not take account of the structural 0 values for within-type relationships,
#  assuming that each element in the bipartite matrix could be related to all N+M-1 other nodes,
#  not just to the other-mode nodes


### Display the first 5 values ---------------------------------------------

SouthernWomen_DegreeWomen[1:5, ]
SouthernWomen_DegreeEvents[1:5, ]


### Most central women and events ---------------------------------------------

#  Let's find the most central women in terms of degree centrality:
#  I will use the head() function to display only the first few elements:

head(sort(SouthernWomen_DegreeWomen[, 1], decreasing = TRUE))

#  And the most central events:

head(sort(SouthernWomen_DegreeEvents[, 1], decreasing = TRUE))


## Betweenness centrality ---------------------------------------------

#  We can compute the betweenness centrality by passing the bipartite network to xBetweennessCentrality().
#  The raw scores will be fine. Normalized ones are problematic 
#  because the standard calculation of the maximum possible betweenness 
#  does not take account of the structural 0 values for within-type relationships.

SouthernWomen_Betweenness <- xBetweennessCentrality(SouthernWomen_Bipartite)

#  Recall from Section 4: 
#  the xUCINET documentation indicates that xBetweennessCentrality() double-counts geodesics 
#  in undirected networks, treating them as if they were directed paths. 
#  This yields undirected counts that are twice as large 
#  as the actual number of geodesics on which a given node lies. 
#  To adjust for this, you can divide the calculated raw scores by 2:

SouthernWomen_Betweenness_adj <- SouthernWomen_Betweenness[, 1] / 2


### Display the first 5 values ---------------------------------------------

#  Note that the object storing the betweenness scores have N + M elements
#  We can separate them between women and events as follows:

SouthernWomen_Betweenness_Women <-
  SouthernWomen_Betweenness_adj[1:NROW(SouthernWomen_Incidence)] # 1 to N elements
SouthernWomen_Betweenness_Events <-
  SouthernWomen_Betweenness_adj[(NROW(SouthernWomen_Incidence) + 1):NROW(SouthernWomen_Bipartite)] # N + 1 to N + M elements

#  Display the first 5 values for each of the vectors above:

SouthernWomen_Betweenness_Women[1:5]
SouthernWomen_Betweenness_Events[1:5]


### Most central women and events ---------------------------------------------

#  Let's find the most central women in terms of betweenness centrality:
#  I will use the head() function to display only the first few elements:

head(sort(SouthernWomen_Betweenness_Women, decreasing = TRUE))

#  And the most central events:

head(sort(SouthernWomen_Betweenness_Events, decreasing = TRUE))




#  4. Visualizing the two-mode network ---------------------------------------------

## Using the incidence matrix ---------------------------------------------

par(mar = c(0, 0, 0, 0))
set.seed(02138)
SouthernWomen_PlotIncidence <- gplot(
  SouthernWomen_Incidence,
  gmode = "twomode",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80"
)

#  With the incidence matrix as the input data and gmode set to "twomode",
#  the function automatically distinguishes the two modes using red circles and blue squares.


## Using the bipartite network ---------------------------------------------

#  Alternatively, we could have used the bipartite network with gmode set to "graph".

set.seed(02138)
SouthernWomen_PlotBipartite <- gplot(
  SouthernWomen_Bipartite,
  gmode = "graph",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80"
)

#  This produces the same node layout as above (provided that the same seed is supplied):

identical(SouthernWomen_PlotIncidence, SouthernWomen_PlotBipartite)

#  This is because, as the function documentation explains, 
#  when gmode is set to "twomode", the supplied two-mode data is converted to a bipartite network
#  before the plotting proceeds. 

#  The disadvantage of using the bipartite network with gmode set "graph"
#  is that the function won't automatically distinguish the modes


## Highlighting the most central nodes ---------------------------------------------

#  We will create a new version of the graph that highlights 
#  the 2 women and the 2 events most central in terms of betweenness
#  Let's first extract the labels for these 2 women and 2 events:

SouthernWomen_MostCentralBetwenness_Women <-
  names(sort(SouthernWomen_Betweenness_Women, decreasing = TRUE)[1:2])

SouthernWomen_MostCentralBetwenness_Events <-
  names(sort(SouthernWomen_Betweenness_Events, decreasing = TRUE)[1:2])

#  We will use different colors to highlight the most central women and events
#  Central women will be colored green, and central events will be colored blue.
#  All other nodes will be colored red

#  We will proceed in two steps. First, we will color all the nodes red
#  by repeating the character "red" as many times as the number of rows in SouthernWomen_Bipartite
#  (the number of women plus the number of events):

SouthernWomen_color <- rep("red", NROW(SouthernWomen_Bipartite))

#  Next, we will replace the colors for the most central women using their row indices in the bipartite network
#  The match() function can be handy here.
#  It returns a vector of the positions of the matches of the first argument in the second argument:

SouthernWomen_MostCentralBetwenness_WomenIndex <-
  match(SouthernWomen_MostCentralBetwenness_Women,
        rownames(SouthernWomen_Bipartite))

SouthernWomen_MostCentralBetwenness_WomenIndex

#  NORA is on the 14th row of SouthernWomen_Bipartite and EVELYN on the 1st

#  Change the color of the 1st and 14th elements of SouthernWomen_color to green:

SouthernWomen_color[SouthernWomen_MostCentralBetwenness_WomenIndex] <- "green"

#  Repeat the steps above for the most central events to color them blue:

SouthernWomen_MostCentralBetwenness_EventsIndex <-
  match(SouthernWomen_MostCentralBetwenness_Events,
        rownames(SouthernWomen_Bipartite))

SouthernWomen_color[SouthernWomen_MostCentralBetwenness_EventsIndex] <- "blue"

#  We can now pass the SouthernWomen_color vector to the vertex.col argument of gplot()
#  (Note that the function will still automatically distinguish the modes using vertex shape,
#  so we don't need to worry about that)

set.seed(02138)
gplot(
  SouthernWomen_Incidence,
  gmode = "twomode",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80",
  vertex.col = SouthernWomen_color,
)




#  5. Core/Periphery partitioning ---------------------------------------------

## Using the xDualCorePeriphery() function ---------------------------------------------

#  We can use the xDualCorePeriphery() function to conduct a core/periphery mapping.
#  It separates the entities in both the row and column sets into two groups such that
#  the core set of row entities is closely linked with the core set of column entities,
#  and the peripheral set of row entities is linked to the core set of column entities, 
#  but not to the peripheral set of column entities

#  To accomplish this, xDualCorePeriphery() passes the one-mode projections to the xCorePeriphery() function.
#  This function, in turn, uses the eigenvector centrality score to measure node coreness. 
#  So it first computes the eigenvector centralities. 
#  Then it finds the best core-periphery split. 
#  It does so by considering all N-1 possible partitions, with core sizes varying from 1 to N-1. 
#  For each possible partition, it computes the correlation between the eigenvector centrality scores and an ideal set of scores 
#  in which the core members have scores of 1 and the peripheral members scores of 0. 
#  Finally, it chooses the partition with the highest correlation

#  The only argument in xDualCorePeriphery() is the incidence matrix:

SouthernWomen_CP <- xDualCorePeriphery(SouthernWomen_Incidence)

#  The function produces an R list containing two vectors, one for the mapping of row entities
#  and the other for the mapping of column entities in the two-mode network.
#  In these, 0s assign entities to the peripheral block while 1s assign them to the core block

#  The plots produced display the normalized eigenvector centrality scores,
#  which are used to produce the core/periphery mapping,
#  arranged in decreasing order for both the row and column entities 
#  (that's why we have two plots)


## Checking membership ---------------------------------------------

#  Let's check the membership of women in the core/periphery:

split(rownames(SouthernWomen_Incidence), SouthernWomen_CP[[1]])

#  Similarly for events:

split(colnames(SouthernWomen_Incidence), SouthernWomen_CP[[2]])


## Subgroup density table ---------------------------------------------

#  We can use the vectors stored in SouthernWomen_CP to create a subgroup density table:

xDensity(SouthernWomen_Incidence,
         ROWS = SouthernWomen_CP[[1]],
         COLS = SouthernWomen_CP[[2]])

#  We can see that the core/periphery partition is moderately successful

#  We can also visualize this in a graph using the plot.mat() function from the blockmodeling package:

library(blockmodeling)
plot.mat(SouthernWomen_Incidence,
         clu = list(SouthernWomen_CP[[1]],
                    SouthernWomen_CP[[2]]))


## Visualizing the core/periphery mapping in a graph ---------------------------------------------

#  Using a procedure similar to the one described above,
#  I will create a graph coloring the core women green and the core events blue
#  The peripheral women and events will be colored red

SouthernWomen_color2 <- rep("red", NROW(SouthernWomen_Bipartite))
SouthernWomen_color2[match(names(which(SouthernWomen_CP[[1]] == 1)), rownames(SouthernWomen_Bipartite))] <- "green"
SouthernWomen_color2[match(names(which(SouthernWomen_CP[[2]] == 1)), rownames(SouthernWomen_Bipartite))] <- "blue"

set.seed(02138)
gplot(
  SouthernWomen_Incidence,
  gmode = "twomode",
  displaylabels = TRUE,
  usearrows = FALSE,
  label.cex = 0.8,
  label.pos = 3,
  edge.col = "grey80",
  vertex.col = SouthernWomen_color2,
)


