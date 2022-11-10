#' Provides the outdegree for multiple types of network relations among the same nodes and a measure of heterogeneity (IQV)
#'
#' Takes an object of class 'list', with a set of different types of network relations (matrices) among the same nodes stored as elements.
#' For the different types of networks it provides for each node:
#'  - the number of valid alters, i.e., the number of alters for which information is available (not equal to NA),
#'  - the degree, i.e., the sum of ties (or tie strengths) to alters (Measures="SumStrength"),
#'  - the average number of ties (or tie strengths) to alters relative to the number of valid cases (Measures="AvStrength"),
#'  - the proportion of ties (or tie strengths) for a specific relation relative to the total ties over all relations (Measures="PropStrength"),
#'  - a measure of heterogeneity using Agresti's index of qualitative variation (IQV) (Measures="IQVType").
#'
#' @param LIST1 A set of matrices stored as a list object
#' @param Measures The measures to be calculated about the set of networks. Options are: "SumStrength","AvStrength","PropStrength" and "IQVType".
#' @param Loops Logical, whether the diagonal values should be ignored (default setting is FALSE, i.e., to ignore the diagonal values)
#' @param Mode Whether the set of networks are one-mode "OneMode" or two-mode networks "TwoMode".
#'
#' @return A matrix with different measures (columns) for each node (row)
#' @export
#' @references Chapter 8 in Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#'
#' @examples
#' ## Consider the different types of networks in the Hawthorne bank wiring room network.
#' ## These are elements 4 to 9 in the project/dataset "Hawthorne_BankWiring".
#' ## All networks are binary, except the last, 9th network ("TradeJobs").
#' ## Although the function can deal with valued networks, we decide to dichotomize the last
#' ## network first and store this as a new (10th) element in the list object:
#' Hawthorne_BankWiring$TradeJobsDich<-(Hawthorne_BankWiring$TradeJobs>0)*1
#' ## We can now apply the function to all binary networks:
#' xMultipleTieComposition(Hawthorne_BankWiring[c(4:8,10)])
#' ## See Chapter 8 in Borgatti, Everett, Johnson and Agneessens (2022)

xMultipleTieComposition<-function(LIST1,
                                  Measures=c("SumStrength","AvStrength","PropStrength","IQVType"),
                                  Loops=FALSE,
                                  Mode="OneMode")
{
  N_NET1<-length(LIST1)
  N_MATRIXk_R<-sapply(LIST1,nrow)
  N_MATRIXk_C<-sapply(LIST1,ncol)
  #==== 2.1. Set up output ====
  LIST2<-LIST1
  for (MATi in 1:N_NET1)
  {
    if (Loops==FALSE) {diag(LIST1[[MATi]])<-NA}
    LIST2[[MATi]]<-!is.na(LIST1[[MATi]])
  }
  OUTPUT1<-sapply(LIST2,rowSums)
  if(is.null(names(LIST1)))
  {
    NAMESMATRICES<-paste("M",c(1:N_NET1),sep="")
  }
  else
  {
    NAMESMATRICES<-names(LIST1)
  }
  LABELS1<-paste("Valid.",NAMESMATRICES,sep="")
  #==== 2.2. Optional measures ====
  if ("SumStrength" %in% Measures)
  {
    OUTPUT2<-sapply(LIST1,rowSums,na.rm=TRUE)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS2<-paste("Sum.",NAMESMATRICES,sep="")
    LABELS1<-c(LABELS1,LABELS2)
  }
  if ("AvStrength" %in% Measures)
  {
    OUTPUT2<-sapply(LIST1,rowMeans,na.rm=TRUE)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS2<-paste("Av.",NAMESMATRICES,sep="")
    LABELS1<-c(LABELS1,LABELS2)
  }
  if ("PropStrength" %in% Measures)
  {
    OUTPUT2<-sapply(LIST1,rowSums,na.rm=TRUE)
    OUTPUT2<-OUTPUT2/apply(as.matrix(OUTPUT2),1,sum)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS2<-paste("Prop.",NAMESMATRICES,sep="")
    LABELS1<-c(LABELS1,LABELS2)
  }
  if ("IQVType" %in% Measures)
  {
    OUTPUT2<-sapply(LIST1,rowSums, na.rm=TRUE)
    OUTPUT2<-OUTPUT2/apply(as.matrix(OUTPUT2),1,sum)
    OUTPUT3<-(1-apply(as.matrix(OUTPUT2^2),1,sum))/(1-1/(ncol(OUTPUT2)))
    OUTPUT1<-cbind(OUTPUT1,OUTPUT3)
    LABELS1<-c(LABELS1,"IQVType")
  }
  #==== 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  return(OUTPUT1)
}


#' Provides different measures about the outgoing ties for nodes on a valued network
#'
#' Takes an object of class 'matrix', which is a valued network, and calculates different measures:
#'  - the sum of the number of valid alters, i.e., the number of alters for which information is available (not equal to NA),
#'  - the sum of tie strengths to alters (Measures="SumStrength"),
#'  - the average tie strengths to alters relative to the number of valid cases (Measures="AvStrength"),
#'  - the standard deviation of tie strengths to alters (Measures="SDStrength"),
#'  - the median of tie strengths to alters (Measures="MedianStrength"),
#'  - the range of tie strengths to alters, i.e., the maximum minus the minimum value (Measures="RangeStrength"),
#'  - the minimum tie strengths to alters (Measures="MinStrength"),
#'  - the maximum tie strengths to alters (Measures="MaxStrength"),
#'  - the interquartile for tie strengths to alters (Measures="IQRStrength"),
#'  - the first quartile for tie strengths to alters (Measures="Q1Strength"),
#'  - the third quartile for tie strengths to alters (Measures="Q3Strength"),
#'  - the IQV for tie strengths to alters (Measures="IQVStrength"),
#'  - the distribution of tie strengths to alters (Measures="DistrStrength").
#'
#' @param NET1 A matrix stored as an object of class matrix
#' @param Measures The measures to be calculated about the set of networks. Options are: "SumStrength","AvStrength","SDStrength","MedianStrength","RangeStrength","MinStrength","MaxStrength","IQRStrength","Q1Strength","Q3Strength","IQVStrength", and "DistrStrength".
#' @param Loops Logical, whether the diagonal values should be ignored (default setting is FALSE, i.e., to ignore the diagonal values)
#' @param Mode Whether the networks is one-mode "OneMode" or two-mode networks "TwoMode".
#' @param Exclude0 Whether to exclude values 0 from the analysis (default is FALSE).
#' @param NrCategoriesIQV An option value to be used for IQV strength, in case there are some additional values/categories that should be considered, but might not be present in the data. A value 1 means that the observed number of categories should be used.
#'
#' @return A matrix with different measures (columns) for each node (row)
#' @export
#' @references Chapter 8 in Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' @importFrom stats IQR median quantile sd
#' @examples
#' ## Applied to the valued network for Zachary's Karate club:
#' xValuedTieComposition(Zachary_KarateClub$Strength, Exclude0=TRUE, NrCategoriesIQV = 20)
#' ## See Chapter 8 in Borgatti, Everett, Johnson and Agneessens (2022)

xValuedTieComposition<-function(NET1,
                                Measures=c("SumStrength","AvStrength","SDStrength",
                                           "MedianStrength","RangeStrength","MinStrength","MaxStrength",
                                           "IQRStrength","Q1Strength","Q3Strength",
                                           "IQVStrength","DistrStrength"),
                                Loops=FALSE,
                                Mode="OneMode",
                                Exclude0=FALSE,
                                NrCategoriesIQV=1)
{
  if(is.matrix(NET1)==F) {stop('Object needs to be of class "matrix"')}
  if(Loops==FALSE) {diag(NET1)<-NA}
  if(Exclude0==TRUE) {NET1[NET1==0]<-NA}

  #==== 2.1. Set up output ====
  #valid cases by node
  OUTPUT1<-matrix(apply(!is.na(NET1),1,sum),nrow(NET1),1)
  LABELS1<-c("Valid")
  #==== 2.2. Optional measures ====
  #Interval variables
  if ("SumStrength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,sum,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"SumStrength")
  }
  if ("AvStrength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,mean,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"AvStrength")
  }
  if ("SDStrength" %in% Measures)
  {
    #    OUTPUT2<-sqrt((OUTPUT1[,1]-1)/OUTPUT1[,1]) * apply(NET1,1,sd,na.rm = T)
    OUTPUT2<-apply(NET1,1,sd,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"SDStrength")
  }
  #Ordinal variables
  if ("MedianStrength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,median,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"MedianStrength")
  }
  if ("RangeStrength" %in% Measures)
  {
    RANGEV<-apply(NET1,1,range,na.rm = T)
    OUTPUT2<-RANGEV[2,]-RANGEV[1,]
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"RangeStrength")
  }
  if ("MinStrength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,min,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"MinStrength")
  }
  if ("MaxStrength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,max,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"MaxStrength")
  }
  if ("IQRStrength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,IQR,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"IQRStrength")
  }
  if ("Q1Strength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,quantile,.25,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Q1Strength")
  }
  if ("Q3Strength" %in% Measures)
  {
    OUTPUT2<-apply(NET1,1,quantile,.75,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Q3Strength")
  }
  if ("IQVStrength" %in% Measures)
  {
    #Idenfity range of integer values
    MaxV<-max(NET1,na.rm = T)
    MinV<-min(NET1,na.rm = T)
    OUTPUT3<-matrix(0,nrow(NET1),MaxV-MinV+1)
    for (cati in MinV:MaxV)
    {
      OUTPUT3[,(cati-MinV+1)]<-apply(NET1==cati,1,sum,na.rm = T)
    }
    OUTPUT3P<-OUTPUT3/apply(OUTPUT3,1,sum)
    if (NrCategoriesIQV==1)
    {
      NCAT<-NCOL(OUTPUT3P)
    } else {
      NCAT<-max(NCOL(OUTPUT3P),NrCategoriesIQV)
      if ((NrCategoriesIQV>1)*(NrCategoriesIQV<NCOL(OUTPUT3P)))
      {
        cat("WARNING: The observed number of unique categories is higher than specified in NrCategoriesIQV", "\n", "\n")
      }
    }
    OUTPUT2<-(1-apply(as.matrix(OUTPUT3P^2),1,sum))/(1-1/(NCAT))
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"IQVStrength")
  }
  if ("DistrStrength" %in% Measures)
  {
    #Idenfity range of integer values
    MaxV<-max(NET1,na.rm = T)
    MinV<-min(NET1,na.rm = T)
    OUTPUT2<-matrix(0,nrow(NET1),MaxV-MinV+1)
    for (cati in MinV:MaxV)
    {
      OUTPUT2[,(cati-MinV+1)]<-apply(NET1==cati,1,sum,na.rm = T)
    }
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS2<-paste("F.",c(MinV:MaxV),sep="")
    LABELS1<-c(LABELS1,LABELS2)
  }
  #==== 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  #restore warnings (1.2)
  return(OUTPUT1)
}

#' Information based on ego's alters values on a categorical variable.
#'
#' This function provide for each ego information about the (relative) presence of alters
#' with specific values on an categorical (nominal) variable and calculates measures of
#' heterogeneity based on these values.
#'
#' The input is a one-mode or two-mode network (\code{NET1}) and a vector containing
#' a categorical (nominal) variable with a (limited) set of values (\code{ATTR1}).
#' It extracts for each node (ego) the nominal values for the alters that
#' ego is directed connected with, and based on this information, provides the
#' frequencies with which the alters of ego belong to each of the (k)
#' unique categorical value.
#' Based on the frequency it then calculates a measure of heterogeinity (IQV - Index of
#' Qualitative Variation).
#'
#' The output contain a table with the frequencies for each categorical value and a
#' measure of heterogeneity for the values of a categorical variable among each ego's alters.
#'
#' As an output it provides by default:
#'  - "Valid": the number of cases for which a valid outgoing connections ego has,
#'  - "EgoValue": what the nominal value for ego is (only possible for one-mode network)
#'  - "Degree": what the "binary" outdegree is for ego, i.e. how many unique outgoing connections ego has: "Degree"
#'  - "W.Degr": the sum of the outgoing ties weighted by their strength.
#'     (For binary network data this is the same as "Degree").
#'  - "F.1", "F.2", ... k columns, where each columns gives the number of times ego is connected with alters with a specific value ("1","2",...)
#'  - "P.1", "P.2", ... k columns, where each columns gives the proportion of ego's alters that belong to category "1", "2", ...
#'  - "BlauH": Blau's measure of heterogeneity, is 1 minus the sum of the squares of the proportions of each value of the categorical variable in ego's network.
#'  A person connected to equal numbers of sociologists, psychologists and economists, will have a Heterogeneity measure of 0.333,  calculated as 1 - ( (1/3)^2 + (1/3)^2) + (1/3)^2) ).
#' \deqn{ BlauH = 1-\sum(p[k]^2)}
#'  - "IQV": Index of Qualitative Variation is a normalized version of this index and is equal to the previous column divided by 1-1/k.
#' Using:
#'  \code{\link{xCreateProject}}: function in this package.
#'
#' @param NET1 A network dataset stored as an object of class "matrix".
#'        A one-mode network with (n) nodes, or a two-mode network
#'        with (m) nodes of mode "A" (rows) and (n) nodes of mode "B" (columns).
#' @param ATTR1 A categorical (nominal) attribute vector of size (n), which
#'        contains (k) unique categories.
#' @param Measures The choice of output to be provided. Choices are:
#'        \itemize{
#'          \item "Count", which provides as many columns as there are unique
#'           categories (k) in the nominal variable (\code{ATTR1}), and gives the number of times
#'           ego is connected with alters with a specific value using ("F.1","F.2",...).
#'          \item "Prop", which provides as many columns as there are unique
#'           categories (k) in the nominal variable (\code{ATTR1}), and gives the proportion of
#'           times ego is connected with alters with a specific value using ("P.1","P.2",
#'           ...).
#'          \item "BlauH", Blau's measure of heterogeneity, is 1 minus the sum of the
#'          squares of the proportions of each value of the categorical variable in ego's
#'          network.
#'          \item "IQV", provides the Index of Qualitative Variation, which is a normalized
#'          version of the Blau index, which is the value for the Blau index divided by 1-1/k.
#'          }
#' @param Loops Whether to ignore the Loops (default is to ignore Loops).
#' @param Mode Whether NET1 is a one-mode ("OneMode") or two-mode ("TwoMode") network.
#'        The default is "Autodetect", which will detect whether the network is a one-mode
#'        or a two-mode network based on the information provided, or based on additional
#'        information linked to the dataset.
#' @return A matrix containing different measures (columns) for each node (row).
#' @export
#'
#' @examples
#' ## Using Krackhardt's friendship network among high-tech managers with department as
#' ## categorical attribute.
#' xAlterCompositionCat(Krackhardt_HighTech$Friendship, Krackhardt_HighTech$Attributes$Department)

xAlterCompositionCat<-function(NET1, ATTR1,
                               Measures=c("Count", "Prop", "BlauH","IQV"),
                               Loops=FALSE,
                               Mode="OneMode")
{
  if(is.matrix(NET1)==FALSE) {stop('Matrix object needs to be of class "matrix"')}
  if(is.vector(ATTR1)==FALSE)  {stop('Attribute object needs to be a vector')}
  if(class(ATTR1)!="numeric" & class(ATTR1)!="integer") {stop('Attribute object needs to contain numeric values')}

  #==== 1.3. Checks data input and options ====
  N_ATTR1<-length(ATTR1)
  N_NET1_R<-nrow(NET1)
  N_NET1_C<-ncol(NET1)

  #= 1.1. Check that length of attribute is equal to the number of columns in matrix:
  if (N_ATTR1!=N_NET1_C)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'ERROR: The number of elements in the attribute vector is ',N_ATTR1, ',', '\n',
         '   and this is different from the number of columns in the matrix, which is ', N_NET1_C, '\n', sep="")
    stop('== == == == Function aborted. Please see error above. == == == == == == == == == == == ==')
  }

  #= 1.2. Check mode and ExcludeD
  #== F. Checks if matrix is two-mode:
  if (Mode=="TwoMode")
  {
    if (N_NET1_R==N_NET1_C)
    {
      cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
           'WARNING: The network was indicated as being two-mode, yet the number of rows and columns are the same.', '\n',
           '   In what follows we will deal with the network as a two-mode network, although this might be a mistake.', '\n')
    }
    if (Loops==FALSE)
    {
      cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
           'WARNING: The network was indicated as being two-mode, yet it was requested that the Loops would be removed.', '\n',
           '   This is not possible for a two-mode network. The Loops was not removed.', '\n')
      Loops<-TRUE
    }
  }
  #== G. Checks if matrix is indicated as one-mode:
  if ((Mode=="OneMode")*(N_NET1_R!=N_NET1_C))
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'WARNING: The network was indicated as being one-mode, yet the number of rows and columns are not the same.', '\n',
         '   In what follows we will deal with the network as a two-mode network.', '\n')
    Mode<-"TwoMode"
  }
  #== B. Check that row and column labels of the matrix correspond:
  if (Mode=="OneMode")
  {
    if (sum(rownames(NET1)!=colnames(NET1))>0)
    {
      cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
           'WARNING: The names in the rows and in the columns of the matrix do not seem to be the same.', '\n',
           '   In what follows we will nevertheless deal with the network as a one-mode network.', '\n')
    }
    if (Loops==FALSE)
    {
      # Remove any diagonal values from the matrix
      diag(NET1)<-NA
    }
  }

  # 2. OUTPUT
  #= 2.1. Set up output ====
  # Create a matrix with alter's attribute values in each column
  AlterCAT<-matrix(ATTR1,NROW(NET1),NCOL(NET1),byrow=T)
  #valid cases by node and node's attribute
  OUTPUT1<-matrix(apply(!is.na(NET1*AlterCAT),1,sum),nrow(NET1),1)
  DEGW<-matrix(apply(NET1,1,sum,na.rm=T),nrow(NET1),1)
  DEG1<-matrix(apply(NET1>0,1,sum,na.rm=T),nrow(NET1),1)

  if (Mode=="OneMode")
  {
    #Add the number of valid cases and attribute vector to the output
    OUTPUT1<-cbind(OUTPUT1,ATTR1,DEG1,DEGW)
    LABELS1<-c("Valid","EgoValue","Degree","W.Degr")
  }

  if (Mode=="TwoMode")
  {
    #Add the number of valid cases and attribute vector to the output
    OUTPUT1<-cbind(OUTPUT1,DEG1,DEGW)
    LABELS1<-c("Valid","Degree","W.Degr")
  }

  #= 2.2. Basic calculations ====
  UNIQ<-unique(ATTR1)
  N_CAT<-length(UNIQ)
  UNIQ<-UNIQ[order(UNIQ)]
  OUTPUT2<-matrix(c(9999),NROW(NET1),N_CAT)
  OUTPUT2
  q1<-1
  for (CAT1 in UNIQ)
  {
    OUTPUT2[,q1]<-apply(NET1*(AlterCAT==CAT1),1,sum,na.rm=T)
    q1<-q1+1
  }

  #= 2.3. Option: Frequency ====
  if ("Count" %in% Measures)
  {
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS2<-paste("F.",UNIQ,sep="")
    LABELS1<-c(LABELS1,LABELS2)
    OUTPUT1
  }
  #= 2.4. Option: Percentage ====
  if ("Prop" %in% Measures)
  {
    OUTPUT3<-OUTPUT2/apply(as.matrix(OUTPUT2),1,sum)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT3)
    LABELS3<-paste("P.",UNIQ,sep="")
    LABELS1<-c(LABELS1,LABELS3)
  }
  #= 2.5. Blau's measure of heterogeneity ====
  if ("BlauH" %in% Measures)
  {
    OUTPUT3<-OUTPUT2/apply(as.matrix(OUTPUT2),1,sum)
    OUTPUT4<-1-apply(as.matrix(OUTPUT3^2),1,sum)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT4)
    LABELS1<-c(LABELS1,"BlauH")
  }
  #= 2.6. Option: IQV ====
  if ("IQV" %in% Measures)
  {
    OUTPUT3<-OUTPUT2/apply(as.matrix(OUTPUT2),1,sum)
    OUTPUT4<-(1-apply(as.matrix(OUTPUT3^2),1,sum))/(1-1/N_CAT)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT4)
    LABELS1<-c(LABELS1,"IQV")
  }
  #= 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  rownames(OUTPUT1)<-rownames(NET1)
  return(OUTPUT1)
}

#' Information based on ego's alters values on a continuous variable.
#'
#' This function takes a one-mode or two-mode matrix (\code{NET1}) and a vector containing a continuous nodal attribute (\code{ATTR1}).
#' It returns for each node in the rows of the matrix information related to the attribute data of this node's alters.
#' For valued networks it weights the impact of alters by the strength of its ties.
#' For a two-mode network, results are provided for the row nodes based on attribute data of the alters (column nodes).
#' In this case the attribute needs to be a vector containing information for the column nodes.
#' As an output it provides by default:
#'  - "SumAlterV": the (weighted) sum of the values for the alters of ego (Sj(Xij*Rj)),
#'  - "AvAlterV": the (weighted) sum of the values for the alters of ego and divides by the number of (non-missing) alters ((Sj(Xij*Rj))/(N-1)),
#'  - "WAvAlterV": the (weighted) average of the values for the alters of ego as a proportion of the strength of the ties ((Sj(Xij*Rj))/(Sj(Xij))),
#'  - "MinAlterV": the minimum value for the alters of ego,
#'  - "MaxAlterV": the maximum value for the alters of ego,
#'  - "RangeAlterV": the range among the values for the alters of ego,
#'  - "SDAlterV": the standard deviation of the values for the alters of ego,
#'
#' @param NET1 A network dataset stored as an object of class "matrix".
#'        A one-mode network with (n) nodes, or a two-mode network
#'        with (m) nodes of mode "A" (rows) and (n) nodes of mode "B" (columns).
#' @param ATTR1 A continuous attribute vector of size (n).
#' @param Measures The choice of output to be provided. See description.
#' @param Loops Whether to ignore the Loops (default is to ignore the Loops).
#' @param Exclude0 Whether to exclude values 0 from the analysis (default is FALSE).
#' @param Mode Whether NET1 is a one-mode ("OneMode") or two-mode ("TwoMode") network.
#' @return A matrix containing different measures (columns) for each node (row).
#' @export
#'
#' @examples
#' ## Consider the Friendship network for Krackhardt's high-tech managers and the attribute tenure.
#' xAlterCompositionCon(Krackhardt_HighTech$Friendship, Krackhardt_HighTech$Attributes$Tenure,
#'                      Exclude0=TRUE)

xAlterCompositionCon<-function(NET1, ATTR1,
                               Measures=c("SumAlterV","AvAlterV","WAvAlterV",
                                          "MinAlterV","MaxAlterV","RangeAlterV","SDAlterV"),
                               Loops=FALSE, Exclude0=FALSE,
                               Mode="OneMode")
{
  if(is.matrix(NET1)==F) {stop('Object needs to be of class "matrix"')}
  if(Loops==FALSE) {diag(NET1)<-NA}
  if(Exclude0==TRUE) {NET1[NET1==0]<-NA}

  N_ATTR1<-length(ATTR1)
  N_NET1_R<-nrow(NET1)
  N_NET1_C<-ncol(NET1)
  #==== 2.1. Set up output ====
  #valid cases by node and node's attribute
  OUTPUT1<-matrix(apply(!is.na(NET1),1,sum),nrow(NET1),1)
  LABELS1<-c("Valid")
  OUTPUT1<-cbind(OUTPUT1,ATTR1)
  LABELS1<-c(LABELS1,"EgoValue")
  AlterValues<-matrix(ATTR1,NROW(NET1),NCOL(NET1),byrow=T)
  AlterValuesC<-AlterValues*NET1
  #==== 2.2. Optional measures ====
  if ("SumAlterV" %in% Measures)
  {
    OUTPUT2<-apply(AlterValuesC,1,sum,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"SumAlterV")
  }
  if ("AvAlterV" %in% Measures)
  {
    OUTPUT2<-apply(AlterValuesC,1,sum,na.rm=T)/apply(!is.na(AlterValuesC),1,sum,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"AvAlterV")
  }
  if ("WAvAlterV" %in% Measures)
  {
    OUTPUT2<-apply(AlterValuesC,1,sum,na.rm=T)/apply(NET1*(!is.na(AlterValuesC)),1,sum,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"WAvAlterV")
  }
  if ("MinAlterV" %in% Measures)
  {
    AlterValuesC0<-AlterValuesC
    AlterValuesC0[NET1==0]<-NA
    OUTPUT2<-apply(AlterValuesC0,1,min,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"MinAlterV")
  }
  if ("MaxAlterV" %in% Measures)
  {
    AlterValuesC0<-AlterValuesC
    AlterValuesC0[NET1==0]<-NA
    OUTPUT2<-apply(AlterValuesC0,1,max,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"MaxAlterV")
  }
  if ("RangeAlterV" %in% Measures)
  {
    AlterValuesC0<-AlterValuesC
    AlterValuesC0[NET1==0]<-NA
    OUTPUT2<-apply(AlterValuesC0,1,max,na.rm=T)-apply(AlterValuesC0,1,min,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"RangeAlterV")
  }
  if ("SDAlterV" %in% Measures)
  {
    WMEAN<-(apply(AlterValuesC,1,sum,na.rm=T)/apply(NET1*(!is.na(AlterValuesC)),1,sum,na.rm=T))
    OUTPUT2<-sqrt((OUTPUT1[,1]-1)/OUTPUT1[,1]) * apply(AlterValuesC,1,sd,na.rm = T)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"SDAlterV")
  }
  #==== 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  rownames(OUTPUT1)<-rownames(NET1)
  return(OUTPUT1)
}


#' Information based the similarity/difference on a categorical between ego and the alters that ego is directly connected with.
#'
#' This function takes a one-mode network (\code{NET1}) and a vector containing a categorical nodal attribute (\code{ATTR1}).
#' It returns for each node in the rows of the matrix information related to the similarity or difference in attribute data between the node and its alters.
#' As an output it provides by default:
#'  -"a": the number of connected alters that belong to the same category as ego.
#'  -"b": the number of connected alters that belong to a different category than ego.
#'  -"c": the number of alters that ego is not connected to, that belong to the same category as ego.
#'  -"b": the number of alters that ego is not connected to, that belong to a different category than ego.
#   -"PctSame": Percentage of alters (weighted by tie strength) that belong to the same type as ego (a/(a+b)).
#   -"EIIndex": Classic measure that takes the number of alters (weighted by tie strength) that belong to a different
#      type than ego (External = b), and subtracts the number of alters that belong to the same type as ego (Internal = a).
#      This difference is divided by the total number of alters (a+b).
#   -"OddsRatio": Takes the odds of alters that belong to the same category (a) divided by alters that belong to a different category (b),
#      compared to the same odds, but for people not connected to ego (c and d).
#   -"LogOddsRatio": Takes the natural log of the OddsRatio
#   -"YulesQ": [(a*d)-(b*c)]/[(a*d)+(b*c)].
#'
#' @param NET1 A network dataset stored as an object of class "matrix".
#'        A one-mode network with (n) nodes, or a two-mode network
#'        with (m) nodes of mode "A" (rows) and (n) nodes of mode "B" (columns).
#' @param ATTR1 A continuous attribute vector of size (n).
#' @param Measures The choice of output to be provided. See description.
#' @param Loops Whether to ignore the Loops (default is to ignore the Loops).
#' @return A matrix containing different measures (columns) for each node (row).
#' @export
#'
#' @examples
#' ## Consider the Friendship network for Krackhardt's high-tech managers and the attribute age.
#' xEgoAlterSimilarityCat(Krackhardt_HighTech$Friendship,Krackhardt_HighTech$Attributes$Department)


xEgoAlterSimilarityCat<-function(NET1, ATTR1,
                                 Measures=c("PctSame","EIIndex","OddsRatio","LogOddsRatio","YulesQ"),
                                 Loops=TRUE)
{
  #==== 1.1. Check all the measures to be included from the full list below ====
  MeasuresAll<-c("PctSame","EIIndex","OddsRatio","LogOddsRatio","YulesQ")
  MeasuresSel<-MeasuresAll %in% Measures
  #Identify any measures that are not among the list of options
  WRONGMEASURES<-(Measures %in% MeasuresAll)==0
  #Give a warning for any unidentified measures
  if (sum(WRONGMEASURES)>0)
  {
    cat("WARNING: The following measure(s) could not be identified, and might contain a typo:", "\n",
        Measures[WRONGMEASURES], "\n", "\n")
  }
  #==== 1.2. Suppress default warnings in output ====
  defaultW <- getOption("warn")
  options(warn = -1)
  #==== 1.3. Checks data input and options ====
  N_ATTR1<-length(ATTR1)
  N_NET1_R<-nrow(NET1)
  N_NET1_C<-ncol(NET1)
  #== A. Check that length of attribute is equal to the number of columns in matrix:
  if (N_ATTR1!=N_NET1_C)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'ERROR: The number of elements in the attribute vector is ',N_ATTR1, ',', '\n',
         '   and this is different from the number of columns in the matrix, which is ', N_NET1_C, '\n', sep="")
    stop('== == == == Function aborted. Please see error above. == == == == == == == == == == == ==')
  }
  #== C. Check that labels of the matrix correspond to those of the attribute file:
  if (sum(rownames(ATTR1)!=colnames(NET1))>0)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'WARNING: The names in the rows of the attribute file do not correspond to those in the matrix.', '\n')
  }
  #== Y. Check if matrix is indeed a one-mode network:
  if (N_NET1_R!=N_NET1_C)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'ERROR: The number of rows is ', N_NET1_R, ', while the number of columns is ', N_NET1_C, '\n',
         '   This function requires a square matrix', '\n', sep="")
    stop('== == == == Function aborted. Please see error above. == == == == == == == == == == == ==')
  }
  diag(NET1)<-NA

  #==== 2.1. Set up output ====
  #valid cases by node and node's attribute
  OUTPUT1<-matrix(apply(!is.na(NET1),1,sum),nrow(NET1),1)
  OUTPUT1<-cbind(OUTPUT1,ATTR1)
  AlterValues<-matrix(ATTR1,nrow(NET1),nrow(NET1),byrow=T)
  SameValues<-(AlterValues==t(AlterValues))
  OtherValues<-(AlterValues!=t(AlterValues))
  SameValuesNC<-apply((NET1==0)*SameValues,1,sum,na.rm=T)
  SameValuesC<-apply((NET1>0)*SameValues,1,sum,na.rm=T)
  OtherValuesNC<-apply((NET1==0)*OtherValues,1,sum,na.rm=T)
  OtherValuesC<-apply((NET1>0)*OtherValues,1,sum,na.rm=T)
  OUTPUT1<-cbind(OUTPUT1,SameValuesC,OtherValuesC,SameValuesNC,OtherValuesNC)
  LABELS1<-c("Valid","EgoType","TieToSame(a)","TieToDiff(b)","NoTieToSame(c)","NoTieToDiff(d)")
  #==== 2.2. Optional measures ====
  if ("PctSame" %in% Measures)
  {
    OUTPUT2<-SameValuesC/(SameValuesC+OtherValuesC)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"PctSame")
  }
  if ("EIIndex" %in% Measures)
  {
    OUTPUT2<-(OtherValuesC-SameValuesC)/(OtherValuesC+SameValuesC)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"EIIndex")
  }
  #IF VALUED...
  VALUEDNET<-(max(NET1, na.rm=T))>1
  if (VALUEDNET==T)
  {
    WSameValues<-apply(NET1*SameValues,1,sum,na.rm=T)
    WOtherValues<-apply(NET1*OtherValues,1,sum,na.rm=T)
    OUTPUT1<-cbind(OUTPUT1,WSameValues,WOtherValues)
    LABELS1<-c(LABELS1,"W.TieToSame(a*)","W.TieToDiff(b*)")
    if ("PctSame" %in% Measures)
    {
      OUTPUT2<-WSameValues/(WSameValues+WOtherValues)
      OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"W.PctSame")
    }
    if ("EIIndex" %in% Measures)
    {
      OUTPUT2<-(WOtherValues-WSameValues)/(WOtherValues+WSameValues)
      OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"W.EIIndex")
    }
  }
  #BACK TO MEASURES
  if ("OddsRatio" %in% Measures)
  {
    OUTPUT2<-(SameValuesC*OtherValuesNC)/(OtherValuesC*SameValuesNC)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"OddsRatio")
  }
  if ("LogOddsRatio" %in% Measures)
  {
    OUTPUT2<-log((SameValuesC*OtherValuesNC)/(OtherValuesC*SameValuesNC))
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"LogOddsRatio")
  }
  if ("YulesQ" %in% Measures)
  {
    OUTPUT2<-((SameValuesC*OtherValuesNC)-(SameValuesNC*OtherValuesC))/
      ((SameValuesC*OtherValuesNC)+(SameValuesNC*OtherValuesC))
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"YulesQ")
  }
  #==== 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  rownames(OUTPUT1)<-rownames(NET1)
  options(warn = defaultW)
  return(OUTPUT1)
}

#' Information based the similarity in continuous characteristics between ego and the alters that ego is directly connected with.
#'
#' This function takes a one-mode network (\code{NET1}) and a vector containing a continuous nodal attribute (\code{ATTR1}).
#' It returns for each node in the rows of the matrix information related to the similarity or difference in attribute data between the node and its alters.
#' For valued networks it weights the impact of alters by the strength of its ties.
#' As an output it provides by default:
#'   -"AvDiffAv": (Weighted) average of difference between alter's values and ego's own value.
#'   -"AvAbsDiff": (Weighted) average of the absolute difference between alter's values and ego's own value.
#'   -"AvSqDiff": (Weighted) average of the squared differences between alter's values and ego's own value.
#'   -"Corr_AbsDiff": For any ego, correlation between ego's (valued) link to alter and the absolute difference in value between ego and alter. Has n-1 alters as cases.
#'   -"Corr_Zegers": For any ego, correlation between ego's link to alter and the product of values for ego and alter divided by their squared values.
#'   -"Corr_MinOverMax": For any ego, correlation between ego's link to alter and the lowest over the highest value.
#'   -"Corr_Product": For any ego, correlation between ego's link to alter and the product of values for ego and alter.
#'
#' @param NET1 A network dataset stored as an object of class "matrix".
#'        A one-mode network with (n) nodes, or a two-mode network
#'        with (m) nodes of mode "A" (rows) and (n) nodes of mode "B" (columns).
#' @param ATTR1 A continuous attribute vector of size (n).
#' @param Measures The choice of output to be provided. See description.
#' @param Loops Whether to ignore the Loops (default is to ignore the Loops).
#' @return A matrix containing different measures (columns) for each node (row).
#' @importFrom stats cor
#' @export
#'
#' @examples
#' ## Consider the Friendship network for Krackhardt's high-tech managers and the attribute age.
#' xEgoAlterSimilarityCon(Krackhardt_HighTech$Friendship,Krackhardt_HighTech$Attributes$Age)

xEgoAlterSimilarityCon<-function(NET1, ATTR1,
                                 Measures=c("AvDiff","AvAbsDiff","AvSqDiff",
                                            "Corr_AbsDiff","Corr_DiffSq","Corr_Zegers","Corr_MinOverMax","Corr_Product"),
                                 Loops=FALSE)
{
  #==== 1.1. Check all the measures to be included from the full list below ====
  MeasuresAll<-c("DiffAv","AbsDiffAv","SqDiff",
                 "Corr_AbsDiff","Corr_DiffSq","Corr_Zegers","Corr_MinOverMax","Corr_Product")
  MeasuresSel<-MeasuresAll %in% Measures
  #Identify any measures that are not among the list of options
  WRONGMEASURES<-(Measures %in% MeasuresAll)==0
  #Give a warning for any unidentified measures
  if (sum(WRONGMEASURES)>0)
  {
    cat("WARNING: The following measure(s) could not be identified, and might contain a typo:", "\n",
        Measures[WRONGMEASURES], "\n", "\n")
  }
  #==== 1.3. Checks data input and options ====
  N_ATTR1<-length(ATTR1)
  N_NET1_R<-nrow(NET1)
  N_NET1_C<-ncol(NET1)
  #== A. Check that length of attribute is equal to the number of columns in matrix:
  if (N_ATTR1!=N_NET1_C)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'ERROR: The number of elements in the attribute vector is ',N_ATTR1, ',', '\n',
         '   and this is different from the number of columns in the matrix, which is ', N_NET1_C, '\n', sep="")
    stop('== == == == Function aborted. Please see error above. == == == == == == == == == == == ==')
  }
  #== C. Check that labels of the matrix correspond to those of the attribute file:
  if (sum(rownames(ATTR1)!=colnames(NET1))>0)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'WARNING: The names in the rows of the attribute file do not correspond to those in the matrix.', '\n')
  }
  #== Y. Check if matrix is indeed a one-mode network:
  if (N_NET1_R!=N_NET1_C)
  {
    cat( '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         'ERROR: The number of rows is ', N_NET1_R, ', while the number of columns is ', N_NET1_C, '\n',
         '   This function requires a square matrix', '\n', sep="")
    stop('== == == == Function aborted. Please see error above. == == == == == == == == == == == ==')
  }
  if (Loops==FALSE) {diag(NET1)<-NA}

  #==== 2.1. Set up output ====
  OUTPUT1<-matrix(apply(!is.na(NET1),1,sum),nrow(NET1),1)
  OUTPUT1<-cbind(OUTPUT1,ATTR1)
  LABELS1<-c("Valid","EgoValue")
  AlterValues<-matrix(ATTR1,NCOL(NET1),NCOL(NET1),byrow=T)
  NET1PERC<-NET1/rowSums(NET1, na.rm=T)
  NET1PERC[NET1==0]<-NA
  OUTPUT2<-rowSums(NET1>0,na.rm=T)
  OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
  LABELS1<-c(LABELS1,"Degree")

  #==== 2.2. Optional measures ====
  if ("AvDiff" %in% Measures)
  {
    OUTPUT2<-rowSums( (AlterValues-t(AlterValues))*NET1PERC , na.rm=T)/OUTPUT1[,3]
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"AvDiff")
  }
  if ("AvAbsDiff" %in% Measures)
  {
    OUTPUT2<-rowSums( abs(AlterValues-t(AlterValues))*NET1PERC , na.rm=T)/OUTPUT1[,3]
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"AvAbsDiff")
  }
  if ("AvSqDiff" %in% Measures)
  {
    OUTPUT2<-rowSums( (AlterValues-t(AlterValues))^2*NET1PERC , na.rm=T)/OUTPUT1[,3]
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"AvSqDiff")
  }
  if ("Corr_AbsDiff" %in% Measures)
  {
    AbsDiff<-abs(t(AlterValues)-AlterValues)
    diag(AbsDiff)<-NA
    OUTPUT2<-mapply(cor,as.data.frame(AbsDiff),as.data.frame(t(NET1)),use="pairwise.complete.obs")
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Corr_AbsDiff")
  }
  if ("Corr_DiffSq" %in% Measures)
  {
    DiffSq<-(t(AlterValues)-AlterValues)^2
    diag(DiffSq)<-NA
    OUTPUT2<-mapply(cor,as.data.frame(DiffSq),as.data.frame(t(NET1)),use="pairwise.complete.obs")
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Corr_DiffSq")
  }
  if ("Corr_Zegers" %in% Measures)
  {
    Zegers<-t(AlterValues)*AlterValues/(t(AlterValues)^2+AlterValues^2)
    diag(Zegers)<-NA
    OUTPUT2<-mapply(cor,as.data.frame(Zegers),as.data.frame(t(NET1)),use="pairwise.complete.obs")
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Corr_Zegers")
  }
  if ("Corr_MinOverMax" %in% Measures)
  {
    MinOverMax<-pmin(t(AlterValues),AlterValues)/pmax(t(AlterValues),AlterValues)
    diag(MinOverMax)<-NA
    OUTPUT2<-mapply(cor,as.data.frame(MinOverMax),as.data.frame(t(NET1)),use="pairwise.complete.obs")
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Corr_MinOverMax")
  }
  if ("Corr_Product" %in% Measures)
  {
    Prod<-t(AlterValues)*AlterValues
    diag(Prod)<-NA
    OUTPUT2<-mapply(cor,as.data.frame(Prod),as.data.frame(t(NET1)),use="pairwise.complete.obs")
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Corr_Product")
  }
  #==== 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  rownames(OUTPUT1)<-rownames(NET1)
  return(OUTPUT1)
  #return(AlterValuesC)
}

#' Measures to capture structural holes.
#'
#' This function takes a one-mode network (\code{NET1}) and calculates measures of structural equivalence.
#' It returns for each node in the rows of the matrix information about effective size and/or constraint.
#' Both measures also work for valued networks.
#'
#' @param NET1 A one-mode network dataset stored as an object of class "matrix".
#' @param Measures The choice of output to be provided. This can be "EffectiveSize", and/or "Constraint".
#' @param Loops Whether to include the Loops (default is to ignore the Loops).
#' @param Include Whether to consider the whole network ("WholeNetwork") or only the egonetwork ("EgoNetwork") around the node.
#' @param Direction How to define the alters for a node (at the moment only outgoing choices).
#' @return A matrix containing different measures (columns) for each node (row).
#' @export
#'
#' @examples
#' ## Consider the Friendship network for Krackhardt's high-tech managers.
#' xStructuralHoles(Padgett_FlorentineFamilies$Marriage)
#' ## Consider the Friendship network for Krackhardt's high-tech managers.
#' xStructuralHoles(Krackhardt_HighTech$Friendship)

xStructuralHoles<-function(NET1,
                           Measures=c("EffectiveSize","Constraint"),
                           Loops=FALSE,
                           Include="WholeNetwork",
                           Direction="Out")
{
  if (Loops==FALSE) {diag(NET1)<-0}
  #==== 2.1. Set up output ====
  OUTPUT1<-matrix(apply(!is.na(NET1),1,sum),nrow(NET1),1)
  LABELS1<-c("Valid")
  OUTPUT2<-matrix(apply(NET1>0,1,sum, na.rm=TRUE),nrow(NET1),1)
  OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
  LABELS1<-c(LABELS1,"Size")

  NET1DEGREE<-rowSums(NET1,na.rm=TRUE)
  NET1P<-NET1/(NET1DEGREE+(NET1DEGREE==0))

  #==== 2.2. Optional measures ====
  if ("EffectiveSize" %in% Measures)
  {
    NET1MAX<-apply(NET1,1,max,na.rm=TRUE)
    NET1M<-NET1/(NET1MAX+(NET1MAX==0))
    DYADICREDUNDANCY<-NET1P%*%t(NET1M)
    diag(DYADICREDUNDANCY)<-0
    OUTPUT2<-rowSums((1-DYADICREDUNDANCY)*(NET1>0), na.rm=TRUE)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"EffectiveSize")
  }

  if ("Constraint" %in% Measures & Include=="WholeNetwork")
  {
    NET1D<-NET1>0
    DYADCONSTR<-((NET1P+(NET1P%*%NET1P))^2)*NET1D
    OUTPUT2<-rowSums(DYADCONSTR)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"ConstraintWhole")
  }
  if ("Constraint" %in% Measures & Include=="EgoNetwork")
  {
    NET1D<-NET1>0
    NET1V<-(NET1D%*%t(NET1))
    DYADCONSTR<-((NET1P+((NET1P/(NET1V+(NET1V==0)))))^2)*NET1D
    OUTPUT2<-rowSums(DYADCONSTR)
    OUTPUT1<-cbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"ConstraintEgo")
  }
  #==== 3.1. Return output ====
  colnames(OUTPUT1)<-LABELS1
  rownames(OUTPUT1)<-rownames(NET1)
  #    return(OUTPUT1)
  # DYADCONSTR
  OUTPUT1
}

