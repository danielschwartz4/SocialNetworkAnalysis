install.packages("haven")
install.packages("devtools")
remove.packages("dplyr")
install.packages("dplyr")
library(xUCINET)
library(sna)
library(haven)
library(tidyverse)
# library(devtools)
# devtools::install_github("tidyverse/dplyr")
library(dplyr)

network <- Corporation$NetworkInfo
attributes <- Corporation$Attributes
corperation <- Corporation$Corporation

directed <- xSymmetrize(corperation, Type = "Max")

# 1. 
EGO <- 3
Alters_Ego3 <- directed[EGO,]>0
Alters_Ego3[EGO] <- 1
Ego3Network <- directed[Alters_Ego3==1,Alters_Ego3==1]
Ego3Network
gplot(Ego3Network, 
      displaylabels = TRUE)

# 2.
Alters_Ego3[EGO] <- 0
Ego3Network <- directed[Alters_Ego3==1,Alters_Ego3==1]
NROW(Ego3Network)
xDensity(Ego3Network)

# 3. 
directedHoles <- xStructuralHoles(directed)
undirectedHoles <- xStructuralHoles(corperation)
directedHoles

# 4. 
# Size
tapply(holes[,2], attributes$seniority, mean)
# EffectiveSize
tapply(holes[,3], attributes$seniority, mean)
# ConstraintWhole
tapply(holes[,4], attributes$seniority, mean)

# 5. 
GSS_network <- read_dta("~/documents/college/SocialNetworkAnalysis/HW5/GSS1985_NetworkExtract.dta", encoding = NULL)

mean(GSS_network$educ1, na.rm = TRUE)
mean(GSS_network$educ2, na.rm = TRUE)
mean(GSS_network$educ3, na.rm = TRUE)
mean(GSS_network$educ4, na.rm = TRUE)
mean(GSS_network$educ5, na.rm = TRUE)
sd(GSS_network$educ1, na.rm = TRUE)
sd(GSS_network$educ2, na.rm = TRUE)
sd(GSS_network$educ3, na.rm = TRUE)
sd(GSS_network$educ4, na.rm = TRUE)
sd(GSS_network$educ5, na.rm = TRUE)

library(tidyverse)
GSS_network <- GSS_network %>%
  mutate(across(
    close12:close45,
    ~ case_when(. == 1 ~ 1,
                . == 2 ~ 0.5,
                . == 3 ~ 0,
                TRUE ~ NA_real_)
  ))
packageVersion("dplyr")

numGivenDir <- xSymmetrize(matrix(GSS_network$numgiven))
numGivenDir
xStructuralHoles(matrix(GSS_network$numgiven))

















