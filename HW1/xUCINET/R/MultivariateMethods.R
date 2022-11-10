
#' Performs multidimensional scaling on a matrix
#'
#' Takes a distance matrix and performs multidimensional scaling on the matrix. This function relies on cmdscale function.
#'
#' @param MAT1 A distance matrix
#' @param Margins The margins added to the output figure
#' @param Distance The distance of the names (labels) of the nodes to the midpoint (as a % of the overall distance).
#'
#' @return A figure containing the MDS solution and the output from cmdscale.
#'
#' @importFrom stats cmdscale
#' @importFrom graphics par points text
#'
#' @export
#'
#' @references
#' Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @examples
#' CITIES1<-matrix(c(
#'      0,  206,  429, 1504,  963, 2976,
#'    206,    0,  233, 1308,  802, 2815,
#'    429,  233,    0, 1075,  671, 2684,
#'   1504, 1308, 1075,    0, 1329, 3273,
#'    963,  802,  671, 1329,    0, 2013,
#'   2976, 2815, 2684, 3273, 2013,    0),6,6)
#' colnames(CITIES1)<-c("BOSTON","NY","DC","MIAMI","CHICAGO","SEATTLE")
#' xMDS(CITIES1,Distance=c(0,4))

xMDS<-function(MAT1, Margins=0.05, Distance=c(0,4))
{
  COOR<-cmdscale(MAT1, eig=TRUE, k = 2)
  par(mar = c(2,2,2,2))
  EXTREME1<-matrix(0,2,2)
  #extra margin for adding names
  MARG<-Margins
  Distance<-Distance/100
  EXTREME1[1,1]<-max(COOR$points[,1])+MARG*sum(abs(range(COOR$points[,1])))
  EXTREME1[2,1]<-min(COOR$points[,1])-MARG*sum(abs(range(COOR$points[,1])))
  EXTREME1[1,2]<-max(COOR$points[,2])+MARG*sum(abs(range(COOR$points[,2])))
  EXTREME1[2,2]<-min(COOR$points[,2])-MARG*sum(abs(range(COOR$points[,2])))
  plot(EXTREME1, xlab=NA, ylab=NA,pch=0,col="white")
  points(COOR$points, xlab=NA, ylab=NA,pch=16,col="grey40")
  text(COOR$points[,1]+Distance[1]*sum(abs(range(COOR$points[,1]))),COOR$points[,2]+Distance[2]*sum(abs(range(COOR$points[,2]))),
       labels=colnames(MAT1))
  COOR
}

#--------------------------------------------------------------------
#' Provides the hierarchical clustering solution for a square matrix
#'
#' Using a square matrix, representing either similarities or differences, this function calculates the clustering solution.
#' This function relies on hclust and plot, and therefore incorporates their respective options. If requested (by adding values to NOC),
#' cluster membership is offered for a solution with k clusters.
#'
#' @param NET1 A binary, undirected one-mode network stored as a 'matrix' object.
#' @param Input Whether the input matrix contains similarities ("Similarities") or differences ("Differences")
#' @param Method which clustering method to use. Default is "complete". See the hclust function.
#' @param TitleDendrogram The title to be used in the plot
#' @param NOC If containing a single value or vector, cluster membership is provided as output.
#'
#' @return A matrix where each row and column are unique cliques found, and the values are the number of nodes two cliques have in common.
#'
#' @importFrom stats as.dist hclust cutree
#'
#' @export
#'
#' @references
#' Chapter 11. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xCliquesMembership()], [xUCINET::xCliquesMembership()], [hclust()]
#'
#' @examples
#' ## Examples of undirected networks (See chapter 11 of Borgatti et al., 2022):
#' CCM<-xCliquesCoMembership(ASNR_Fig11x1,Min=3)
#' xHierarchicalClustering(CCM, Input="Similarities", Method="single")
#' MATY<-matrix(c(2,3,1,
#'                3,5,6,
#'                1,6,4),3,3)
#' xHierarchicalClustering(MATY, Input="Similarities")

xHierarchicalClustering<-function(NET1, Input=NULL, Method="complete",
                                  TitleDendrogram = "Cluster Dendrogram", NOC=NULL)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Matrix NET1 needs to be of class "matrix" ###')}
  if(is.null(Input)) {stop(' .  ### Please specify whether the input NET1 contains similarities or differences (Input="Similarities") or (Input="Differences") ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(NR1!=NC1){stop(' .   ### Matrix NET1 needs to be symmetrical. ')}

  if(Input=="Similarities"){
    DISTM<-stats::as.dist(-NET1)
    YLAB<-"High <-  Similarities (reversed)  -> Low"
    NOTE1<-"NOTE: Y axis represent the negative value of the original SE values."
  }
  if(Input=="Differences"){
    DISTM<-stats::as.dist(NET1)
    YLAB<-"Low <-  Differences  -> High"
    NOTE1<-""
  }
  HCLUST<-stats::hclust(DISTM, method=Method)
  plot(HCLUST, ylab = YLAB,
       main=TitleDendrogram, xlab=deparse(substitute(MAT1)),
       sub = paste(NOTE1),
       cex.axis=1.1, cex.main=1, cex.lab=1.1,
       font.axis=1, col.sub = "red")

  if(!is.null(NOC))
  {
    OUTPUT1<-matrix(NA,NR1,length(NOC))
    for (k in 1:length(NOC))
    {
      OUTPUT1[,k]<-stats::cutree(HCLUST,NOC[k])
    }
    rownames(OUTPUT1)<-rownames(NET1)
    colnames(OUTPUT1)<-paste("CL_",NOC, sep="")
    OUTPUT1
  }
}

#' Performs a correspondence analysis (CA) on a matrix
#'
#' Takes a rectangular matrix and performs a correspondence analysis on the matrix. This function relies on the corresp() function (MASS).
#'
#' @param MAT1 A rectangular matrix
#' @param DIM The number of dimensions to be considered
#' @param Margins The margins added to the output figure
#' @param Distance The distance of the names (labels) of the nodes to the midpoint (as a % of the overall distance).
#' @param Color A vector of size 2, indicating the colors for row and column variables.
#'
#' @return A figure containing the correspondence analysis and the output from corresp().
#'
#' @importFrom MASS corresp
#' @importFrom graphics par points text
#'
#' @export
#'
#' @references
#' Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @examples
#' MATY<-matrix(c(6,3,3,2,
#'                8,2,5,2,
#'                4,8,4,5,
#'                9,8,6,3,
#'                7,3,2,5,
#'                3,2,1,5,
#'                8,5,6,7,
#'                3,2,4,2),8,4,byrow=TRUE)
#' rownames(MATY)<-c("P1","P2","P3","P4","P5","P6","P7","P8")
#' colnames(MATY)<-c("V1","V2","V3","V4")
#' xCorrespondenceAnalysis(MATY)

xCorrespondenceAnalysis<-function(MAT1, Margins=0.05, Distance=c(0,4), DIM=2, Color=c("Blue","DarkRed"))
{
  COOR<-MASS::corresp(MAT1, nf=DIM)
  par(mar = c(2,2,2,2))
  EXTREME1<-matrix(0,2,2)
  #extra margin for adding names
  MARG<-Margins
  COOR2<-rbind(COOR$rscore,COOR$cscore)
  EXTREME1[1,1]<-max(COOR2[,1])+MARG*sum(abs(range(COOR2[,1])))
  EXTREME1[2,1]<-min(COOR2[,1])-MARG*sum(abs(range(COOR2[,1])))
  EXTREME1[1,2]<-max(COOR2[,2])+MARG*sum(abs(range(COOR2[,2])))
  EXTREME1[2,2]<-min(COOR2[,2])-MARG*sum(abs(range(COOR2[,2])))
  plot(EXTREME1, xlab=NA, ylab=NA,pch=0,col="white")
  points(COOR$rscore, xlab=NA, ylab=NA,pch=16,col=Color[1])
  text(COOR$rscore[,1]+Distance[1]/100*sum(abs(range(COOR$rscore[,1]))),COOR$rscore[,2]+Distance[2]/100*sum(abs(range(COOR$rscore[,2]))),
       labels=rownames(MAT1), col=Color[1])
  points(COOR$cscore, xlab=NA, ylab=NA,pch=17,col=Color[2])
  text(COOR$cscore[,1]+Distance[1]/100*sum(abs(range(COOR$cscore[,1]))),COOR$cscore[,2]+Distance[2]/100*sum(abs(range(COOR$cscore[,2]))),
       labels=colnames(MAT1), col=Color[2])
  COOR
}


