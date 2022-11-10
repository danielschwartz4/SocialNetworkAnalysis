#' Dichotomizes the ties of a network
#'
#' This function offers different options for dichotomizing a network. The default is "GT" with Value=0.
#'
#' @param NET1 A network to be dichotomized.
#' @param Value The value at which to dichotomize. When using (Type="EQ") this can also be a vector of values to be set to 1.
#' @param Type The type of dichotomization. The five options are: "GT" (greater than), "LT" (less than), "EQ" (equal),
#   "GTE" (greater than or equal), "LTE" (less than or equal).
#'
#' @return A dichotomized network.
#' @export
#'
#' @examples
#' M1<-matrix(c(1,2,3,4,5,6,7,8,NA),3,3)
#' xDichotomize(M1,Value=c(4,7,2,20),Type="EQ")
#' xDichotomize(M1,Value=c(4),Type="GT")

xDichotomize<-function(NET1, Value=0, Type="GT")
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]

  #Five options: "GT" (greater than), "LT" (less than), "EQ" (equal),
  #"GTE" (greater than or equal), "LTE" (less than or equal)
  if(Type=="GT") {
    OUTPUT1<-(NET1>Value)*1
  }

  if(Type=="LT") {
    OUTPUT1<-(NET1<Value)*1
  }

  if(Type=="EQ") {
    OUTPUT1<-matrix((NET1 %in% Value)*1,NR1,NC1)
    OUTPUT1[is.na(NET1)]<-NA
    rownames(OUTPUT1)<-rownames(NET1)
    colnames(OUTPUT1)<-colnames(NET1)
  }

  if(Type=="GTE") {
    OUTPUT1<-(NET1>=Value)*1
  }

  if(Type=="LTE") {
    OUTPUT1<-(NET1<=Value)*1
  }

  OUTPUT1
}

#' Transposes the ties of a network
#'
#' This function transposes a network containing network data.
#'
#' @param NET1 A network to be transposed.
#'
#' @return The transposed network.
#' @export
#'
#' @examples
#' M1<-matrix(c(1,2,3,4,5,6,7,8,NA),3,3)
#' xTranspose(M1)

xTranspose<-function(NET1)
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]

  #Transpose the network
  OUTPUT1<-t(NET1)
  OUTPUT1
}

#' Symmetrizes the ties of a network
#'
#' This function allows different forms of symmetrization on a network.
#'
#' @param NET1 A network to be symmetrized.

#' @param Type The type of symmetrization. The five options are: "Min" (take the minimum of both ties), "Max" (take the maximum of both ties),
#'  "Av" (take the average of both ties), "Sum" (take the sum of both ties), "Prod" (take the product of both ties).
#'
#' @return The symmetrized network.
#' @export
#'
#' @examples
#' M1<-matrix(c(1,2,3,4,5,6,7,8,NA),3,3)
#' xSymmetrize(M1)

xSymmetrize<-function(NET1, Type="Min")
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]

  #Symmetrize the network

  if(Type=="Min") {
    OUTPUT1<-pmin(NET1,t(NET1),na.rm=TRUE)
    TRANSF<-" minimum of values in both directions"
  }

  if(Type=="Max") {
    OUTPUT1<-pmax(NET1,t(NET1),na.rm=TRUE)
    TRANSF<-" maximum of values in both directions"
  }

  if(Type=="Av") {
    OUTPUT1<-(NET1+t(NET1))/2
    TRANSF<-" average of values in both directions"
  }

  if(Type=="Sum") {
    OUTPUT1<-(NET1+t(NET1))
    TRANSF<-" sum of values in both directions"
  }

  if(Type=="Prod") {
    OUTPUT1<-(NET1*t(NET1))
    TRANSF<-" product of values in both directions"
  }
  OUTPUT1
}

#' Combine two networks
#'
#' This function allows different forms of combining two networks.
#'
#' @param NET1 The first network to be used.
#' @param NET2 The second network to be used.
#' @param Type The way to combine two matrices. The six options are: "Min" (take the minimum of both ties), "Max" (take the maximum of both ties - the default),
#'  "Av" (take the average of both ties), "Sum" (take the sum of both ties), "Prod" (take the product of both ties), "Unique" (creates a unique code for both - only to be used for matrices with limited number of values).
#'
#' @return The combined network.
#' @export
#'
#' @examples
#' M1<-matrix(c(1,2,3,4,5,6,7,8,NA),3,3)
#' M2<-matrix(c(3,2,1,4,5,6,7,8,NA),3,3)
#' xCombineTies(M1,M2)

xCombineTies<-function(NET1, NET2, Type="Max")
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  if(!is.matrix(NET2)) {stop(' .  ### Network file NET2 needs to be of class "matrix" ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  NR2<-dim(NET2)[1]
  NC2<-dim(NET2)[2]
  if(NR1!=NR2) {stop(' .  ### Network file NET1 is not the same dimensions as network file NET2 ###')}
  if(NC1!=NC2) {stop(' .  ### Network file NET1 is not the same dimensions as network file NET2 ###')}

  #Combine two matrices

  if(Type=="Min") {
    OUTPUT1<-pmin(NET1,NET2,na.rm=TRUE)
  }

  if(Type=="Max") {
    OUTPUT1<-pmax(NET1,NET2,na.rm=TRUE)
  }

  if(Type=="Av") {
    OUTPUT1<-(NET1+NET2)/2
  }

  if(Type=="Sum") {
    OUTPUT1<-(NET1+NET2)
  }

  if(Type=="Prod") {
    OUTPUT1<-(NET1*NET2)
  }

  if(Type=="Unique") {
    MAX1<-max(NET1)
    OUTPUT1<-(NET1+(MAX1+1)*NET2)
  }
  OUTPUT1
}

#' Impute missing data
#'
#' This function imputes missing data in a networks.
#'
#' @param NET1 The network to be used.
#' @param Type The way to impute missing data. The current option is: "Density" (take the observed density as a basis).
#' @param Loops Whether to include self-nominations.
#'
#' @return An imputed network.
#' @importFrom stats runif
#' @export
#'
#' @examples
#' M1<-matrix(c(0,NA,1,1,0,0,0,NA,0),3,3)
#' xImputeMissingData(M1, Type="Density")

xImputeMissingData<-function(NET1, Type="Density", Loops=FALSE)
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  if(sum(is.na(NET1))==0) {stop(' .  ### Network file NET1 does not contain any missing values (indicated by NA) ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if((NR1!=NC1) & (Loops==FALSE)) {stop(' .  ### Network file NET1 is not one-mode, and therefore putting Loops=FALSE does not make sense ###')}
  #Impute network

  if(Type=="Density") {
    if(Loops==FALSE) {
      NET1x<-NET1
      diag(NET1x)<--1000
      NMiss<-sum(is.na(NET1x))
      diag(NET1x)<-NA
      NET1[is.na(NET1x)]<-(runif(NMiss,0,1)<mean(NET1x,na.rm=TRUE))
    } else {
      NMiss<-sum(is.na(NET1))
      NET1[is.na(NET1)]<-(runif(NMiss,0,1)<mean(NET1,na.rm=TRUE))
    }
    OUTPUT1<-NET1
  }
  OUTPUT1
}

#' Calculates the geodesic distance for a network
#'
#' This function calculates the geodesic distance for a networks.
#'
#' @param NET1 The network to be used.
#'
#' @return The geodesic distance for a network.
#' @export
#'
#' @examples
#' M1<-matrix(c(0,NA,1,1,0,0,0,NA,0),3,3)
#' xGeodesicDistance(M1)

xGeodesicDistance<-function(NET1)
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]

  #Use igraph
  NET1i<-igraph::graph_from_adjacency_matrix(NET1)
  OUTPUT1<-igraph::distances(NET1i, mode="out")
  OUTPUT1
}

#' Converts an attribute into a network
#'
#' This function converts an attribute vector into a network.
#'
#' @param ATT1 The attribute vector to be used.
#' @param Type The conversion to be used. The options are: "Same" (gives 1 for ij if i and j have the same value),
#'  "AbsDiff" (where ij is the absolute difference in value between i and j),
#'  "Diffij" (where ij is the difference: i - j),
#'  "Diffji" (where ij is the difference: j - i),
#'  "Sender" (where ij is the value for i),
#'  "Received" (where ij is the value for j).
#'
#' @return An network.
#' @export
#'
#' @examples
#' A1<-c(1,3,1,2,2,3,1,3)
#' xAttributeToNetwork(A1, Type="Same")

xAttributeToNetwork<-function(ATT1, Type="Same")
{
  options(scipen = 100)
  # ATTR CHECK
  if(!is.vector(ATT1)) {stop(' .  ### Attribute file ATT1 needs to be of class "vector" ###')}

  NN1<-length(ATT1)
  MATT1<-matrix(ATT1,NN1,NN1)

  if(Type=="Same") {
    OUTPUT1<-(MATT1==t(MATT1))*1
  }

  if(Type=="AbsDiff") {
    OUTPUT1<-abs(MATT1-t(MATT1))
  }

  if(Type=="Diffij") {
    OUTPUT1<-MATT1-t(MATT1)
  }

  if(Type=="Diffji") {
    OUTPUT1<-t(MATT1)-MATT1
  }

  if(Type=="Sender") {
    OUTPUT1<-t(MATT1)
  }

  if(Type=="Receiver") {
    OUTPUT1<-MATT1
  }
  OUTPUT1
}

#' Normalizes the values in a network
#'
#' This function allows different forms of normalizing a network.
#'
#' @param NET1 A network to be normalized.
#'
#' @param Type The type of symmetrization. The five options are:
#'  "Max" (divides the value of each cell value by the maximum of the tie values in a row),
#'  "Av" (divides the value of each cell value by the average of the tie values in a row),
#'  "Sum" (divides the value of each cell by the sum of the tie values in a row),
#'  "Center" (subtracts the mean for a row from each cell value).
#'
#' @return The normalized network.
#' @export
#'
#' @examples
#' M1<-matrix(c(1,2,3,4,5,6,7,8,9),3,3)
#' xNormalize(M1, Type="Max")

xNormalize<-function(NET1, Type)
{
  options(scipen = 100)
  # NETWORK MATRIX CHECK
  if(!is.matrix(NET1)) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}

  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]

  #Symmetrize the network

  if(Type=="Max") {
    OUTPUT1<-NET1/apply(NET1, FUN=max, MARGIN=1, na.rm=TRUE)
    TRANSF<-" divides by maximum of values in a row"
  }

  if(Type=="Av") {
    OUTPUT1<-NET1/rowMeans(NET1,na.rm=TRUE)
    TRANSF<-" divides by average of values in a row"
  }

  if(Type=="Sum") {
    OUTPUT1<-NET1/rowSums(NET1,na.rm=TRUE)
    TRANSF<-" divides by sum of values in a row"
  }

  if(Type=="Center") {
    OUTPUT1<-NET1-rowMeans(NET1,na.rm=TRUE)
    TRANSF<-" subtracts mean in a row from each value"
  }
  OUTPUT1
}
