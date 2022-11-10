library(data.table)
library(xUCINET)

network <- fread("~/documents/college/SocialNetworkAnalysis/hw2/Hansell_network.csv")
attributes <- fread("~/documents/college/SocialNetworkAnalysis/hw2/Hansell_attrib.csv")
# network <- read.csv(file = "~/documents/college/SocialNetworkAnalysis/hw2/Hansell_network.csv", header = TRUE)
# attributes <- read.csv("file = ~/documents/college/SocialNetworkAnalysis/hw2/Hansell_attrib.csv")


project <- xCreateProject(
  GeneralDescription = "Digraph for Student Example",
  NetworkName = "StudentGraph",
  NETFILE1= "~/documents/college/SocialNetworkAnalysis/hw2/Hansell_network.csv",
  FileType = "csv",
  NetworkDescription = "Binary Connections among strong components", 
  Mode = c("Novice"),
  Directed = TRUE,
  Loops = FALSE,
  Values = "Binary",
  References="No references")


# In degree
xDegreeCentrality(project$StudentGraph)

# in degree distribution
dist <- table(xDegreeCentrality(project$StudentGraph)[,1])
barplot(dist,  ylim = c(0,6), main = "In degree distribution", ylab = "Count", xlab = "In Degree")

# Reciprocity
xReciprocity(project$StudentGraph)

project_attr <- data.frame(Gender = c(rep("Boy", times = 12), rep("Girl", times = 15)))
rownames(project_attr) <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27")

project <- xAddAttributesToProject(
  ProjectName = project,
  ATTFILE1 = "~/documents/college/SocialNetworkAnalysis/hw2/Hansell_attrib.csv",
  FileType = "csv",
  AttributesDescription = c("Gender"),
  Mode = c("Novice")
)

project$Attributes
xDensity(project$StudentGraph, ROWS = project$Attributes$Gender, COLS = project$Attributes$Gender)











