library(xUCINET)
library(sna)

network <- Corporation$NetworkInfo
attributes <- Corporation$Attributes
corperation <- Corporation$Corporation


# Subgroup density
xDensity(corperation, ROWS = attributes$seniority, COLS = attributes$seniority)


centrality <- xDegreeCentrality(corperation)
centrality

#  par(mar = c(0, 0, 0, 0))
par(mar= c(1, 1, 1, 1))
set.seed(02138)
gplot(corperation, 
      displaylabels = TRUE, 
      vertex.cex = 0.3 * centrality,
      vertex.col = attributes$seniority + 2, 
      vertex.sides = attributes$office + 2)


fine = 500 # this will adjust the resolving power.
pal = colorRampPalette(c('blue','red'))

#this gives you the colors you want for every point
graphCol = pal(fine)[as.numeric(cut(centrality,breaks = fine))]

# Best design
gplot(corperation, 
      displaylabels = TRUE, 
      vertex.cex = (.13 * attributes$projects)+1,
      vertex.col = graphCol,
      vertex.sides = attributes$seniority+3)
