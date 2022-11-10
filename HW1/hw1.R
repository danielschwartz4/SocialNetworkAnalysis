data <- Zachary_KarateClub$Connection
data

xDensity(data)

gd <- xGeodesicDistance(data)
gd
max(gd)
mean(gd)

xComponents(data)

xTransitivity(data)

xDegreeCentrality(data)
