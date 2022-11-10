#' Calculates the average or sum of all pairs of nodes in a network or for subparts of a network.
#'
#' This function takes the average (Option="Mean") or the sum (Option="Sum") of all (observed) cells in a one-mode or two-mode network.
#' - For a binary (0/1) network taking the sum gives the total number of ties present (value 1) in the network.
#' - For a binary (0/1) network the average is equal to the density, i.e. the proportion of ties with values 1.
#' - For a valued network (a network with weights) taking the sum gives the total weights of all cells in the matrix.
#' - For a valued network taking the average means taking the total of all values in the matrix and divide this by the number of possible ties, which gives the average value for the observed cells in a network. Note that this includes values 0.
#' - Cells with missing data (indicated by values 'NA' in the matrix) are ignored in both the numerator and denominator. Extra information about missing data is provided when running the function.
#' - For one-mode networks you can also indicate whether to ignore the Loops entries (self-nominations). The default is to ignore the Loops (Loops=FALSE).
#' - If an attribute file is specified for rows (ROWS) and/or the columns (COLS), then the average or sum is provided for each unique group in that attribute file. This might be useful to examine whether specific types of actors tend to have specific ties with specific other nodes. Note that ROWS and COLS does not have to contain the same attribute.
#' - If ROWS or COLS is set to "Unique" (default is "All"), then the matrix is not aggregated for that row or column, and so values are obtained for each node. This might be useful if you want to know, for example, if specific nodes have more connections with specific types of actors.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Loops Whether to ignore the Loops (default is to ignore Loops, i.e., Loops=FALSE).
#' @param ROWS Whether to calculate the sum or average over all nodes in the row (="All"), whether to keep each node in the row as a unique entity ("Unique") or whether to aggregate based on some attribute variable (in which case (ROWS= name of the attribute dataset).
#' @param COLS Same as for rows, but here for columns. Note that ROWS and COLS can contain a different attribute variable, for example to see whether there is any difference in male or females bully high achieving versus low achieving students (here gender would be ROWS and academic achievement would be COLS).
#' @param Option Whether to obtain the average ("Mean") or sum ("Sum").
#'
#' @return A single value or matrix calculating the sum or average across all nodes, or split by category (or node).
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xDegreeCentrality()], [xUCINET::xConnectedness()], [xUCINET::xCompactness()], [xUCINET::xReciprocity()]
#'
#' @examples
#'
#' ## A one-mode network
#' xDensity(Krackhardt_HighTech$Advice)
#' xDensity(Krackhardt_HighTech$Advice, ROWS=Krackhardt_HighTech$Attributes$Department,
#'    COLS=Krackhardt_HighTech$Attributes$Department)
#' xDensity(Krackhardt_HighTech$Advice, ROWS="All",
#'    COLS=Krackhardt_HighTech$Attributes$Department, Option="Sum")
#'
#' ## A two-mode network
#' xDensity(Davis_SouthernWomen$Attendance)

xDensity<-function(NET1, Loops=FALSE, ROWS="All", COLS="All", Option="Mean")
{
  options(scipen = 100)
  if(!Option %in% c("Mean", "Sum")) {stop(' .  ### The argument (Option) should be "Mean" or "Sum" ###')}

  #NETWORK MATRIX CHECK
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(Loops==FALSE & NR1==NC1 ) {diag(NET1)<-NA}
  if(Loops==FALSE & NR1!=NC1)
  {
    cat( '\n',' .  == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         ' .   ## WARNING about argument: (Loops=FALSE)', '\n',
         ' .   ## The input is a two-mode network, and therefore the default', '\n',
         ' .   ## setting of removing the "Loops" values is ignored.', '\n',
         ' .  == == == == == == == == == == == == == == == == == == == == == ==', '\n')
  }
  #---------------------------------------------------------------------------------------------------
  # UNIQUE CODES (1=AGGREGATE ALL, -1 KEEP UNIQUE, all others check for attribute file)
  if(length(ROWS)==1)
  {
    if(ROWS=="All") {AG_R<-1}
    if(ROWS=="Unique") {AG_R<-(-1)}
  }
  else {AG_R<-0}

  if(length(COLS)==1)
  {
    if(COLS=="All") {AG_C<-1}
    if(COLS=="Unique") {AG_C<-(-1)}
  }
  else {AG_C<-0}
  #---------------------------------------------------------------------------------------------------
  #SITUATION 1: NO ATTRIBUTES
  if(AG_R==1 & AG_C==1)
  {
    SUMF<-sum(NET1, na.rm=T)
    SUMFVAL<-sum(!is.na(NET1))
  }
  #---------------------------------------------------------------------------------------------------
  #SITUATION 0: both Unique - Returns original
  if(AG_R==-1 & AG_C==-1)
  {
    SUMF<-NET1
    SUMFVAL<-1
  }
  #---------------------------------------------------------------------------------------------------
  #SITUATION 0b: Unique and All
  if(AG_R==-1 & AG_C==1)
  {
    SUMF<-t(t(rowSums(NET1,na.rm=TRUE)))
    SUMFVAL<-t(t(rowSums(!is.na(NET1),na.rm=TRUE)))
    if(is.null(rownames(NET1))) {NAMESR<-c(1:NC1)} else {NAMESR<-rownames(NET1)}
    rownames(SUMF)<-NAMESR
    colnames(SUMF)<-"All"
  }
  #---------------------------------------------------------------------------------------------------
  #SITUATION 0c: Unique and All
  if(AG_R==1 & AG_C==-1)
  {
    SUMF<-t(colSums(NET1,na.rm=TRUE))
    SUMFVAL<-t(colSums(!is.na(NET1),na.rm=TRUE))
    if(is.null(colnames(NET1))) {NAMESC<-c(1:NC1)} else {NAMESC<-colnames(NET1)}
    colnames(SUMF)<-NAMESC
    rownames(SUMF)<-"All"
  }

  #---------------------------------------------------------------------------------------------------
  #SITUATION 2: AT LEAST ONE ATTRIBUTE
  if(AG_R==0)
  {
    #ATTRIBUTES CHECK (IF APPLICABLE)
    if(!is.matrix(ROWS) & !is.vector(ROWS)) {stop('\n',' .  ###  Attribute file for rows (ROWS) needs to be of class "matrix" or "vector"','\n')}
    if(is.matrix(ROWS)){
      if(min(dim(ROWS))!=1) {stop('\n',' .  ###  Attribute file for rows (ROWS) needs to be a single row or column','\n')}
    }
    if(is.matrix(ROWS)) {NAR<-max(dim(ROWS))} else {NAR<-length(ROWS)}
    if(NAR!=NR1) {stop('\n',' .  ###  Attribute file for rows (ROWS) needs to have the same length as the number of rows of the matrix (NET1)','\n')}
  }
  if(AG_C==0)
  {
    #ATTRIBUTES CHECK (IF APPLICABLE)
    if(!is.matrix(COLS) & !is.vector(COLS)) {stop('\n',' .  ###  Attribute file for columns (COLS) needs to be of class "matrix" or "vector"','\n')}
    if(is.matrix(COLS)){
      if(min(dim(COLS))!=1) {stop('\n',' .  ###  Attribute file for columns (COLS) needs to be a single row or column','\n')}
    }
    if(is.matrix(COLS)) {NAC<-max(dim(COLS))} else {NAC<-length(COLS)}
    if(NAC!=NC1) {stop('\n',' .  ###  Attribute file for columns (COLS) needs to have the same length as the number of columns of the matrix (NET1)','\n')}
  }

  #CONSIDER ALL POSSIBILITIES:
  if(AG_R==0 & AG_C==0)
  {
    SUM1<-sapply(by(NET1,ROWS,FUN=colSums, na.rm=T),identity)
    SUMF<-sapply(by(SUM1,COLS,FUN=colSums, na.rm=T),identity)

    SUM1VAL<-sapply(by(!is.na(NET1),ROWS,FUN=colSums, na.rm=T),identity)
    SUMFVAL<-sapply(by(SUM1VAL,COLS,FUN=colSums, na.rm=T),identity)
  }
  if(AG_R==0 & AG_C==-1)
  {
    if(is.null(colnames(NET1))) {NAMESC<-c(1:NC1)} else {NAMESC<-colnames(NET1)}
    SUMF<-t(sapply(by(NET1,ROWS,FUN=colSums, na.rm=T),identity))
    colnames(SUMF)<-NAMESC
    SUMFVAL<-t(sapply(by(!is.na(NET1),ROWS,FUN=colSums, na.rm=T),identity))
  }
  if(AG_R==0 & AG_C==1)
  {
    SUM1<-t(sapply(by(NET1,ROWS,FUN=colSums, na.rm=T),identity))
    SUMF<-t(t(rowSums(SUM1, na.rm=TRUE)))
    SUM1VAL<-t(sapply(by(!is.na(NET1),ROWS,FUN=colSums, na.rm=T),identity))
    SUMFVAL<-t(t(rowSums(SUM1VAL, na.rm=TRUE)))
    colnames(SUMF)<-c("All")
  }
  if(AG_R==-1 & AG_C==0)
  {
    if(is.null(rownames(NET1))) {NAMESR<-c(1:NR1)} else {NAMESR<-rownames(NET1)}
    SUMF<-sapply(by(t(NET1),COLS,FUN=colSums, na.rm=T),identity)
    rownames(SUMF)<-NAMESR
    SUMFVAL<-sapply(by(t(!is.na(NET1)),COLS,FUN=colSums, na.rm=T),identity)
  }
  if(AG_R==1 & AG_C==0)
  {
    SUM1<-sapply(by(t(NET1),COLS,FUN=colSums, na.rm=T),identity)
    SUMF<-t(colSums(SUM1, na.rm=TRUE))
    SUM1VAL<-sapply(by(t(!is.na(NET1)),COLS,FUN=colSums, na.rm=T),identity)
    SUMFVAL<-t(colSums(SUM1VAL, na.rm=TRUE))
    rownames(SUMF)<-c("All")
  }
  if(Option=="Mean") {OUTPUT1<-SUMF/SUMFVAL}
  if(Option=="Sum") {OUTPUT1<-SUMF}

  # A bit of info about missing cells:
  NOTMISSING<-sum(!is.na(NET1))
  if(Loops==FALSE & NR1==NC1) {NPOSS<-NR1*(NR1-1)} else {NPOSS<-NR1*NC1}
  cat( '\n',' .  ------------------------------------------------------------------ ', '\n',
       ' .   Number of valid cells:', NOTMISSING,'\n',
       ' .   which corresponds to:', NOTMISSING/(NPOSS)*100,'% of considered cells.','\n',
       ' .  ------------------------------------------------------------------ ', '\n\n')

  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates the level of reciprocity in a binary, one-mode network.
#'
#' This function takes a one-mode network and provides information related to the level of reciprocity in the network.
#' The basis for calculating the level of reciprocity are the number of mutual (M), asymmetric (A) and null (N) dyads.
#' This is better known as the dyad census or MAN distribution.
#' We can either count the number of dyads connected by a tie (which may or may not be reciprocated)
#' and calculate the proportion of dyads that have reciprocated ties (the dyad based method)
#' or we can count the number of arcs (directed edges) and calculate the proportion of arcs that
#' are reciprocated (the arc based method).
#' Note that the Loops entries (self-nominations) are always ignored for this function and that two-mode data are not allowed.
#' The output generates:
#' - Number and proportion of valid dyads. Dyads where either cell i,j or j,i are missing (indicated by values 'NA' in the matrix) are ignored.
#' - The density for the network, i.e. the proportion of ties with values 1 (considering only those dyads that have a valid value in both directions),
#' - The number of mutual dyads (where (i->j)=1 and (j->i)=1),
#' - The number of asymmetric dyads (where (i->j)=1 and (j->i)=0, or (i->j)=0 and (j->i)=1),
#' - The number of null dyads (where (i->j)=0 and (j->i)=0),
#' - The dyad based reciprocity: M/(M+A)
#' - The arc based reciprocity: 2M/(2M+A), which can be compared to the density
#'
#' @param NET1 A one-mode network stored as a 'matrix' object.#'
#'
#' @return A vector with a series of different measures. See Description.
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xTransitivity()], [xUCINET::xCyclicality()]
#'
#' @examples
#' ## A network with 4 valid dyads, where the density for the values in these 4 dyads is: 5/8.
#' NET1<-matrix(c(0,1,NA,0,
#'       1,0,1,NA,
#'       1,0,1,1,
#'       0,NA,1,0),4,4,byrow=TRUE)
#' xReciprocity(NET1)
#'
#' ## Krackhardt's high-tech managers:
#' xReciprocity(Krackhardt_HighTech$Advice)
#' xReciprocity(Krackhardt_HighTech$Friendship)
#' xReciprocity(Krackhardt_HighTech$ReportTo)
#'
#' ## To show that the arc reciprocity index would be similar to the density for a
#' ## random network, let's generate 100 networks of size 50 by 50 with density .3:
#' ## Store the results for each generated network in RES1:
#' RES1<-matrix(NA,100,8)
#' ## Now generate 100 networks use a loop:
#' for (k1 in 1:100)
#'   {
#'   MA1<-matrix(runif(2500)>=.7,50,50)
#'   RES1[k1,]<-t(xReciprocity(MA1))
#'   }
#' ## We can check that the density for the generated networks is around .3:
#' mean(RES1[,3])
#' ## Now, we can check that the arc reciprocity index for the generated graphs.
#' ## This is indeed around .3:
#' mean(RES1[,8])
#' ## Proving that the arc reciprocity index will be close to the density for a random graph.

xReciprocity<-function(NET1)
{
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ## Object needs to be of class "matrix"')}
  if(NR1!=NC1) {stop(' .  ## Object is not a one-mode network. Number of rows and columns are different.')}
  diag(NET1)<-NA
  if(sum(NET1!=0 & NET1!=1, na.rm=TRUE)>0) {stop(' .  ## Matrix should be binary. The only values allowed are 0, 1 and NA.')}

  OUTPUT1<-matrix(NA,8,1)
  OUTPUT1[1,]<-sum(!is.na(NET1*t(NET1)))/2
  OUTPUT1[2,]<-sum(!is.na(NET1*t(NET1)))/(NR1*(NR1-1))*100
  OUTPUT1[3,]<-mean(NET1*t(NET1<5),na.rm=T)
  OUTPUT1[4,]<-sum(NET1*t(NET1), na.rm=T)/2
  OUTPUT1[5,]<-sum((1-NET1)*t(NET1), na.rm=T)
  OUTPUT1[6,]<-sum((1-NET1)*(1-t(NET1)), na.rm=T)/2
  OUTPUT1[7,]<-OUTPUT1[4,]/(OUTPUT1[4,]+OUTPUT1[5,])
  OUTPUT1[8,]<-OUTPUT1[4,]*2/(OUTPUT1[4,]*2+OUTPUT1[5,])
  colnames(OUTPUT1)<-c("Value")
  rownames(OUTPUT1)<-c("Valid Dyads","%Valid Dyads","Density (valid dyads)","Mutual","Asymmetric","Null","DyadReciprocity","ArcReciprocity")
  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates the level of transitivity in a binary, one-mode network.
#'
#' This function takes a one-mode network and provides information related to the level of transitivity in the network.
#' For an undirected network this is also sometimes referred to as clustering.
#' The basis for calculating the level of transitivity is the number of transitivity triplets (where i->j=1, j->k=1 and i->k=1)
#' divided by the number of opportunities for a transitivity triplet (where i->j=1 and j->k=1, but where i->k could be either 0 or 1).
#' Note that the Loops entries (self-nominations) are always ignored for this function and that two-mode data are not allowed.
#' Triplets where information on any of the three cells (i,j), (j,k) or (i,k) are missing (indicated by values 'NA' in the matrix) are ignored.
#' The output generates:
#' - Number and proportion of valid triplets, i.e., where (i,j), (j,k) and (i,k) are either 0 or 1, but not NA.
#' - The density for the network, i.e. the proportion of ties with values 1 (considering only those triplets that are valid),
#' - The number of triplets belonging to any of the unique 8 different types of triplets,
#' - The Transitivity Index: (ij=1 & jk=1 & ik=1)/((ij=1 & jk=1 & ik=1)+(ij=1 & jk=1 & ik=0)), which can be compared to the density.
#'
#' @param NET1 A one-mode network stored as a 'matrix' object.
#'
#' @return A vector with a series of measures. See Description.
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xReciprocity()], [xUCINET::xCyclicality()]
#'
#' @examples
#' ## A network with 6 valid triplets out of 4*3*2 = 0.25, where the density is 0.722
#' ## and the index is 1/(1+2)=.33.
#' NET1<-matrix(c(0,1,NA,0,
#'       1,0,1,NA,
#'       1,0,1,1,
#'       0,NA,1,0),4,4,byrow=TRUE)
#' xTransitivity(NET1)
#'
#' ## Transitivity for Krackhardt's high-tech managers (index of 0.66, 0.46 and 0.00 respectively):
#' xTransitivity(Krackhardt_HighTech$Advice)
#' xTransitivity(Krackhardt_HighTech$Friendship)
#' xTransitivity(Krackhardt_HighTech$ReportTo)
#'
#' ## To show that the transitivity index would be similar to the density for a
#' ## random network, let's generate 100 networks of size 50 by 50 with density .3:
#' ## Store the results for each generated network in RES1:
#' RES1<-matrix(NA,100,12)
#' ## Now generate 100 networks use a loop:
#' for (k1 in 1:100)
#'   {
#'   MA1<-matrix(runif(2500)>=.7,50,50)
#'   RES1[k1,]<-t(xTransitivity(MA1))
#'   }
#' ## We can check that the density for the generated networks is around .3:
#' mean(RES1[,3])
#' ## Now, we can check that the transitivity index for the generated graphs.
#' ## This is indeed around .3:
#' mean(RES1[,12])
#' ## Proving that the transitivity index will be close to the density for a random graph.

xTransitivity<-function(NET1)
{
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ## Object needs to be of class "matrix"')}
  if(NR1!=NC1) {stop(' .  ## Object is not a one-mode network. Number of rows and columns are different.')}
  diag(NET1)<-NA
  if(sum(NET1!=0 & NET1!=1, na.rm=TRUE)>0) {stop(' .  ## Matrix should be binary. The only values allowed are 0, 1 and NA.')}

  NET1V1<-NET1
  NET1V0<-1-NET1
  NET1V1[is.na(NET1)]<-0
  NET1V0[is.na(NET1)]<-0
  diag(NET1V1)<-0
  diag(NET1V0)<-0
  NET1V<-NET1V1+NET1V0

  OUTPUT1<-matrix(NA,12,1)
  NET1V
  OUTPUT1[1,]<-sum((NET1V%*%NET1V)*NET1V)
  OUTPUT1[2,]<-sum((NET1V%*%NET1V)*NET1V)/(NR1*(NR1-1)*(NR1-2))
  OUTPUT1[4,]<-sum((NET1V1%*%NET1V1)*NET1V1)
  OUTPUT1[5,]<-sum((NET1V1%*%NET1V1)*NET1V0)
  OUTPUT1[6,]<-sum((NET1V1%*%NET1V0)*NET1V1)
  OUTPUT1[7,]<-sum((NET1V0%*%NET1V1)*NET1V1)
  OUTPUT1[8,]<-sum((NET1V0%*%NET1V0)*NET1V1)
  OUTPUT1[9,]<-sum((NET1V1%*%NET1V0)*NET1V0)
  OUTPUT1[10,]<-sum((NET1V0%*%NET1V1)*NET1V0)
  OUTPUT1[11,]<-sum((NET1V0%*%NET1V0)*NET1V0)
  OUTPUT1[12,]<-OUTPUT1[4,]/(OUTPUT1[5,]+OUTPUT1[4,])
  OUTPUT1[3,]<-(OUTPUT1[4,]*3+OUTPUT1[5,]*2+OUTPUT1[6,]*2+OUTPUT1[7,]*2+OUTPUT1[8,]+OUTPUT1[9,]+OUTPUT1[10,])/(3*OUTPUT1[1,])
  rownames(OUTPUT1)<-c("Valid triplets","%Valid triplets","Density (valid triplets)",
                       "ij=1 & jk=1 & ik=1* ",
                       "ij=1 & jk=1 & ik=0* ",
                       "ij=1 & jk=0 & ik=1 ",
                       "ij=0 & jk=1 & ik=1 ",
                       "ij=0 & jk=0 & ik=1 ",
                       "ij=1 & jk=0 & ik=0 ",
                       "ij=0 & jk=1 & ik=0 ",
                       "ij=0 & jk=0 & ik=0 ",
                       "*TransitivityIndex")

  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates the level of cyclicality in a binary, one-mode network.
#'
#' This function takes a one-mode network and provides information related to the level of cyclicality in the network.
#' For an undirected network this is also sometimes referred to as clustering (and will give the same result as transitivity).
#' The basis for calculating the level of cyclicality is the number of cyclical triplets (where i->j=1, j->k=1 and k->i=1)
#' divided by the number of opportunities for a cyclical triplet (where i->j=1 and j->k=1, but where k->i could be either 0 or 1).
#' Note that the Loops entries (self-nominations) are always ignored for this function and that two-mode data are not allowed.
#' Triplets where information on any of the three cells (i,j), (j,k) or (k,i) are missing (indicated by values 'NA' in the matrix) are ignored.
#' The output generates:
#' - Number and proportion of valid triplets, i.e., where (i,j), (j,k) and (k,i) are either 0 or 1, but not NA.
#' - The density for the network, i.e. the proportion of ties with values 1 (considering only those triplets that are valid),
#' - The number of triplets belonging to any of the unique 8 different types of triplets,
#' - The Cyclicality Index: (ij=1 & jk=1 & ki=1)/((ij=1 & jk=1 & ki=1)+(ij=1 & jk=1 & ki=0)), which can be compared to the density.
#'
#' @param NET1 A one-mode network stored as a 'matrix' object.
#'
#' @return A vector with a series of measures. See Description.
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xReciprocity()], [xUCINET::xTransitivity()]
#'
#' @examples
#' ## A network with 6 valid triplets out of 4*3*2 = 0.25, where the density is 0.833
#' ## and the index is 3/(3+1)=.75.
#' NET1<-matrix(c(0,1,NA,0,
#'       1,0,1,NA,
#'       1,0,1,1,
#'       0,NA,1,0),4,4,byrow=TRUE)
#' xCyclicality(NET1)
#'
#' ## Cyclicality for Krackhardt's high-tech managers (index of 0.38, 0.28 and 0.00 respectively):
#' xTransitivity(Krackhardt_HighTech$Advice)
#' xTransitivity(Krackhardt_HighTech$Friendship)
#' xTransitivity(Krackhardt_HighTech$ReportTo)
#'
#' ## To show that the cyclicality index would be similar to the density for a
#' ## random network, let's generate 100 networks of size 50 by 50 with density .3:
#' ## Store the results for each generated network in RES1:
#' RES1<-matrix(NA,100,12)
#' ## Now to generate 100 networks use a loop:
#' for (k1 in 1:100)
#'   {
#'   MA1<-matrix(runif(2500)>=.7,50,50)
#'   RES1[k1,]<-t(xCyclicality(MA1))
#'   }
#' ## We can check that the density for the generated networks is around .3:
#' mean(RES1[,3])
#' ## Now, we can check that the cyclicality index for the generated graphs.
#' ## This is indeed around .3:
#' mean(RES1[,12])
#' ## Proving that the cyclicality index will be close to the density for a random graph.

xCyclicality<-function(NET1)
{
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ## Object needs to be of class "matrix"')}
  if(NR1!=NC1) {stop(' .  ## Object is not a one-mode network. Number of rows and columns are different.')}
  diag(NET1)<-NA
  if(sum(NET1!=0 & NET1!=1, na.rm=TRUE)>0) {stop(' .  ## Matrix should be binary. The only values allowed are 0, 1 and NA.')}

  NET1V1<-NET1
  NET1V0<-1-NET1
  NET1V1[is.na(NET1)]<-0
  NET1V0[is.na(NET1)]<-0
  diag(NET1V1)<-0
  diag(NET1V0)<-0
  NET1V<-NET1V1+NET1V0

  OUTPUT1<-matrix(NA,12,1)
  NET1V
  OUTPUT1[1,]<-sum((NET1V%*%NET1V)*t(NET1V))
  OUTPUT1[2,]<-sum((NET1V%*%NET1V)*t(NET1V))/(NR1*(NR1-1)*(NR1-2))
  OUTPUT1[4,]<-sum((NET1V1%*%NET1V1)*t(NET1V1))
  OUTPUT1[5,]<-sum((NET1V1%*%NET1V1)*t(NET1V0))
  OUTPUT1[6,]<-sum((NET1V1%*%NET1V0)*t(NET1V1))
  OUTPUT1[7,]<-sum((NET1V0%*%NET1V1)*t(NET1V1))
  OUTPUT1[8,]<-sum((NET1V0%*%NET1V0)*t(NET1V1))
  OUTPUT1[9,]<-sum((NET1V1%*%NET1V0)*t(NET1V0))
  OUTPUT1[10,]<-sum((NET1V0%*%NET1V1)*t(NET1V0))
  OUTPUT1[11,]<-sum((NET1V0%*%NET1V0)*t(NET1V0))
  OUTPUT1[12,]<-OUTPUT1[4,]/(OUTPUT1[5,]+OUTPUT1[4,])
  OUTPUT1[3,]<-(OUTPUT1[4,]*3+OUTPUT1[5,]*2+OUTPUT1[6,]*2+OUTPUT1[7,]*2+OUTPUT1[8,]+OUTPUT1[9,]+OUTPUT1[10,])/(3*OUTPUT1[1,])
  rownames(OUTPUT1)<-c("Valid triplets","%Valid triplets","Density (valid triplets)",
                       "ij=1 & jk=1 & ik=1* ",
                       "ij=1 & jk=1 & ik=0* ",
                       "ij=1 & jk=0 & ik=1 ",
                       "ij=0 & jk=1 & ik=1 ",
                       "ij=0 & jk=0 & ik=1 ",
                       "ij=1 & jk=0 & ik=0 ",
                       "ij=0 & jk=1 & ik=0 ",
                       "ij=0 & jk=0 & ik=0 ",
                       "*CyclicalityIndex")

  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates information related to the number of components in a network.
#'
#' Identifies the component for a one-mode or two-mode binary network and provides additional information.
#' It provides information about the number of components, the size of the largest (maximum) and smallest (minimum) component, the mean and median size of the different components.
#' It also calculates the component ratio, which is the largest component minus one divided by n-1.
#' The minimum component size to be considered can be defined by (MinCompSize). The default is to set it to value 1, so isolates are also considered in the component count.
#' For a two-mode network, separate information is provided about both modes.
#' This function relies on the (component.dist) function in the sna package, and for directed networks this offers a choice between ("strong","weak","unilateral", and "recursive").
#' Missing data (indicated by NA) are not allowed.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param Measures Optional measures to be calculated. These include providing information about the size of the smallest component ("Min"), the largest component ("Max"), the component ratio ("ComponentRatio"), the average component size over all components ("Mean"), and the median ("Median").
#' @param MinCompSize The minimum size of the components to be considered. Components smaller than this will be ignored from the analysis. (Default is to consider all, including isolates, i.e., components of size 1).
#' @param Type Whether to consider a weak, strong (or other) way to identify components. This relies on the component.dist function in the sna package, and feeds it directly into the option "connected" in the component.dist function.
#'
#' @importFrom sna component.dist
#' @importFrom stats median
#'
#' @return A matrix containing multiple measures (rows) and a single column (for a one-mode matrix) or 3 columns ("ModeA", "ModeB" and "Both", for a two-mode network).
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#' Carter T. Butts (2020). sna: Tools for Social Network Analysis. R package version 2.6. http://CRAN.R-project.org/package=sna
#'
#' @seealso [sna::component.dist()], [xUCINET::xConnectedness()]
#'
#' ## Examples of undirected one-mode networks (see Borgatti et al., 2022, Figure 10.3)
#' xComponents(ASNR_Fig10x3A)
#' xComponents(ASNR_Fig10x3B)
#' xComponents(ASNR_Fig10x3C)
#' xComponents(ASNR_Fig10x3D)
#' xComponents(ASNR_Fig10x3E)
#'
#' ## An examples of a directed one-mode networks
#' NET1<-matrix(c(0,1,1,0,0,0,0,0,
#'               1,0,1,1,0,0,0,0,
#'               1,1,1,0,0,0,0,0,
#'               0,1,1,0,0,0,0,0,
#'               0,0,0,0,0,0,0,0,
#'               0,0,0,0,0,0,1,0,
#'               0,0,0,0,0,1,0,1,
#'               0,0,0,0,0,0,0,0),8,8,byrow=TRUE)
#' xComponents(NET1, Type="weak")
#' xComponents(NET1, Type="strong")
#' xComponents(NET1, Type="strong",MinCompSize = 2)
#'
#' ## A two-mode network
#' NET2<-matrix(c(0,1,0,0,0,0,0,
#'                1,0,0,0,0,0,0,
#'                1,1,1,0,0,0,0,
#'                0,1,1,0,0,0,0,
#'                0,0,0,1,1,0,0,
#'                0,0,0,0,0,1,0,
#'                0,0,0,0,0,1,0,
#'                0,0,0,0,0,0,0),8,7,byrow=TRUE)
#' xComponents(NET2)
#' xComponents(NET2, MinCompSize = 2)
#' xComponents(NET2, MinCompSize = 3)
#'
xComponents<-function(NET1, Measures=c("Min", "Max", "ComponentRatio", "Mean", "Median"), MinCompSize=1, Type="strong")
{
  #NETWORK MATRIX CHECK
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1))>0) {stop(' .  ## Matrix should not contain NAs')}
  if(sum(NET1!=0 & NET1!=1)>0) {stop(' .  ## Matrix should be binary (0 and 1 values only)')}

  cat( '\n',' .  == == == == == == == == == == == == == == == == == == == ', '\n',
       ' .   Minimum component size considered:', MinCompSize, '\n',
       ' .  == == == == == == == == == == == == == == == == == == == ', '\n\n')
  if(NR1==NC1)
  {
    diag(NET1)<-0
    CSize<-sna::component.dist(NET1,connected=Type)$csize
    CSize<-CSize[CSize>=MinCompSize]
    OUTPUT1<-sum(CSize)
    LABELS1<-"Number of nodes"
    OUTPUT1<-rbind(OUTPUT1,length(CSize))
    LABELS1<-c(LABELS1,"Number of components")
    #==== Optional measures ====
    if ("Min" %in% Measures)
    {
      OUTPUT2<-min(CSize)
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Min componentsize")
    }
    if ("Max" %in% Measures)
    {
      OUTPUT2<-max(CSize)
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Max componentsize")
    }
    if ("ComponentRatio" %in% Measures)
    {
      OUTPUT2<-(length(CSize)-1)/(sum(CSize)-1)
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"ComponentRatio")
    }
    if ("Mean" %in% Measures)
    {
      OUTPUT2<-mean(CSize)
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Mean of componentsizes")
    }
    if ("Median" %in% Measures)
    {
      OUTPUT2<-median(CSize)
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Median of componentsizes")
    }
    colnames(OUTPUT1)<-"Value"
    rownames(OUTPUT1)<-LABELS1
  }

  if(NR1!=NC1)
  {
    #becomes a bipartite automatically unless the number of rows and columns is the same.
    CSizeM<-sna::component.dist(NET1)$membership
    # Remove any components of size smaller than MinCompSize (put to NA)
    CSizeF1<-table(CSizeM)
    GroupsToInclude<-names(CSizeF1[CSizeF1>=MinCompSize])
    IncludeM<-CSizeM %in% (as.numeric(GroupsToInclude))
    CSizeM[IncludeM==0]<-NA

    CSizeF<-table(CSizeM)
    CSizeF_A<-table(CSizeM[1:NR1])
    CSizeF_B<-table(CSizeM[(NR1+1):(NR1+NC1)])

    OUTPUT1<-c(sum(CSizeF_A),sum(CSizeF_B),sum(CSizeF))
    LABELS1<-"Number of nodes for mode A, B and both"

    OUTPUT2<-c(length(CSizeF_A),length(CSizeF_B),length(CSizeF))
    OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
    LABELS1<-c(LABELS1,"Number of components")
    if ("Min" %in% Measures)
    {
      OUTPUT2<-c(min(CSizeF_A),min(CSizeF_B),min(CSizeF))
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Min componentsize")
    }
    if ("Max" %in% Measures)
    {
      OUTPUT2<-c(max(CSizeF_A),max(CSizeF_B),max(CSizeF))
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Max componentsize")
    }
    if ("ComponentRatio" %in% Measures)
    {
      OUTPUT2<-c((max(CSizeF_A)-1)/(OUTPUT1[1,1]-1),(max(CSizeF_B)-1)/(OUTPUT1[1,2]-1),(max(CSizeF)-1)/(OUTPUT1[1,3]-1))
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"ComponentRatio")
    }
    if ("Mean" %in% Measures)
    {
      OUTPUT2<-c(mean(CSizeF_A),mean(CSizeF_B),mean(CSizeF))
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Mean of componentsizes")
    }
    if ("Median" %in% Measures)
    {
      OUTPUT2<-c(median(CSizeF_A),median(CSizeF_B),median(CSizeF))
      OUTPUT1<-rbind(OUTPUT1,OUTPUT2)
      LABELS1<-c(LABELS1,"Median of componentsizes")
    }
    colnames(OUTPUT1)<-c("ModeA","ModeB","Both")
    rownames(OUTPUT1)<-LABELS1
  }
  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates the connectedness (fragmentation) of a network or for subparts of a network.
#'
#' Using a one-mode or two-mode binary network, this function calculates the connectedness (or its inverse the fragmentation of a network).
#' Connectedness is the proportion of pairs of nodes that can reach each other, irrespective of how long it takes, which is equivalent to the proportion of pairs of nodes that are located in the same component.
#' In practice it takes the reachability/connectedness matrix where a cell has value 1 if i and j can reach each other and takes the average (Option="Mean") or the sum (Option="Sum") of pairs of nodes that can reach each other in a network.
#' - For a two-mode network, the result can be split into three parts (within mode A, within mode B, and between mode A and mode B nodes), and therefore we will obtain a matrix with within and between values.
#' - Missing data NA are not allowed.
#' - In addition, if an attribute file is specified for rows (ROWS) and/or the columns (COLS), then the average or sum is provided for each unique group in that attribute file. This might be useful to examine whether specific types of actors tend to be reachable from specific other nodes. Note that ROWS and COLS does not have to contain the same attribute.
#' - If ROWS or COLS is set to "Unique" (default is "All"), then the matrix is not aggregated for that row or column, and so reachability is obtained for each node. This might be useful if you want to know, for example, if specific nodes are reachable from nodes with specific types of actors.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param ROWS Whether to calculate the sum or average of connectedness over all nodes in the row (="All"), whether to keep each node in the row as a unique entity ("Unique") or whether to aggregate based on some attribute variable (in which case (ROWS= name of the attribute dataset).
#' @param COLS Same as for rows, but here for columns. Note that ROWS and COLS can contain a different attribute variable, for example to see whether there is any difference in male or females are able to reach more high achieving versus low achieving students (here gender would be ROWS and academic achievement would be COLS).
#' @param Option Whether to obtain the average ("Mean") or sum ("Sum") for the connectedness matrix.
#'
#' @return A single value or matrix calculating the sum or average of reachability across all nodes, or split by category (or node).
#' @importFrom sna geodist
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' Carter T. Butts (2020). sna: Tools for Social Network Analysis. R package version 2.6. http://CRAN.R-project.org/package=sna
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xCompactness()], [sna::geodist()]
#'
#' @examples
#' ## A simple undirected one-mode network
#' xConnectedness(ASNR_Fig10x3D)
#' ## Adding an attribute for rows and columns:
#' xConnectedness(ASNR_Fig10x3D,ROWS=c(1,1,1,2,2,2),COLS=c(1,1,1,2,2,2))
#'
#' ## Some more complex examples with a two-mode network
#' NET2<-matrix(c(0,1,0,0,0,0,0,
#'               1,0,0,0,0,0,0,
#'               1,1,1,0,0,0,0,
#'               0,1,1,0,0,0,0,
#'               0,0,0,1,1,0,0,
#'               0,0,0,0,0,1,0,
#'               0,0,0,0,0,1,0,
#'               0,0,0,0,0,0,0),8,7,byrow=TRUE)
#' ## When running the function on the twomode, note the use of "A" and "B" to separate
#' ## mode A (row) from mode B (column) nodes:
#' xConnectedness(NET2)
#' ## When adding an attribute for ROWS, note that mode A nodes are now split further
#' ## into 2 groups: A_1, A_2 and A_3.
#' xConnectedness(NET2,ROWS=c(1,1,1,2,2,2,3,3))
#' ## We can also keep the values separate for each unique node in the column:
#' xConnectedness(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS="Unique")
#' ## Or have a different attribute for mode A and mode B:
#' xConnectedness(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS=c("P","E","T","T","P","T","E"))
#' ## Or combine unique (in the rows) with aggregating all (in the columns):
#' xConnectedness(NET2,ROWS="Unique")
#'
xConnectedness<-function(NET1, ROWS="All", COLS="All", Option="Mean")
{
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  if(!Option %in% c("Mean", "Sum")) {stop(' .  ### The argument (Option) should be "Mean" or "Sum" ###')}
  if(sum(is.na(NET1))>0) {stop(' .  ## Matrix should not contain NAs')}
  if(sum(NET1!=0 & NET1!=1)>0) {stop(' .  ## Matrix should be binary (0 and 1 values only)')}

  #NETWORK MATRIX CHECK
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  #CALCULATION
  GDISTNI<-!is.infinite((sna::geodist(NET1))$gdist)

  #SPECIAL WAY TO DEAL WITH TWO-MODE -> BIPARTITE
  if(NR1!=NC1)
  {
    if(length(ROWS)==1)
    {
      if(ROWS=="All") {AG_R<-1}
      if(ROWS=="Unique") {AG_R<-(-1)}
    }
    else {AG_R<-0}

    if(length(COLS)==1)
    {
      if(COLS=="All") {AG_C<-1}
      if(COLS=="Unique") {AG_C<-(-1)}
    }
    else {AG_C<-0}
    # NOW ADD A AND BE BEFORE
    if(AG_R==1){PART_A<-rep("A",NR1)}
    if(AG_R==-1 & is.null(rownames(NET1))) {PART_A<-paste("A",c(1:NR1),sep="_a")}
    if(AG_R==-1 & !is.null(rownames(NET1))) {PART_A<-paste("A",rownames(NET1),sep="_")}
    if(AG_R==0) {PART_A<-paste("A",ROWS,sep="_")}

    if(AG_C==1){PART_B<-rep("B",NC1)}
    if(AG_C==-1 & is.null(colnames(NET1))) {PART_B<-paste("B",c(1:NC1),sep="_e")}
    if(AG_C==-1 & !is.null(colnames(NET1))) {PART_B<-paste("B",colnames(NET1),sep="_")}
    if(AG_C==0) {PART_B<-paste("B",COLS,sep="_")}

    ROWS<-c(PART_A,PART_B)
    COLS<-c(PART_A,PART_B)
  }
  OUTPUT1<-xDensity(GDISTNI, Loops=FALSE, ROWS=ROWS, COLS=COLS, Option=Option)
  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates the compactness (breadth) of a network or for subparts of a network.
#'
#' Using a one-mode or two-mode binary network, this function calculates the compactness (or its inverse the breadth of a network).
#' Compactness is based on the reciprocal of the geodesic distance between pairs of nodes, so that two nodes at distance 1 have a reciprocal value of 1/1, two nodes at distance 2, a value of 1/2, etc.
#' Nodes that cannot reach each other are defined as having a distance infinity, which means the pair will obtain a value 1/inf=0 in the reciprocal distance matrix.
#' In practice it takes this reciprocal distance and takes the average (Option="Mean") or the sum (Option="Sum") of the values for all pairs in the network.
#' - For a two-mode network, the result can be split into three parts (within mode A, within mode B, and between mode A and mode B nodes), and therefore we will obtain a matrix with within and between values.
#' - Missing data NA are not allowed.
#' - In addition, if an attribute file is specified for rows (ROWS) and/or the columns (COLS), then the average or sum is provided for each unique group in that attribute file. This might be useful to examine whether specific types of actors tend on average to be at shorter distances from specific other nodes. Note that ROWS and COLS does not have to contain the same attribute.
#' - If ROWS or COLS is set to "Unique" (default is "All"), then the matrix is not aggregated for that row or column, and so the reciprocal distance to other (sets of) nodes is obtained for each node. This might be useful if you want to know, for example, if specific nodes are at shorter distance from nodes with specific types of actors.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param ROWS Whether to calculate the sum or average of the reciprocal distances over all nodes in the row (="All"), whether to keep each node in the row as a unique entity ("Unique"), or whether to aggregate based on some attribute variable (in which case (ROWS= name of the attribute dataset).
#' @param COLS Same as for rows, but here for columns. Note that ROWS and COLS can contain a different attribute variable, for example to see whether there is any difference in male or female students in their reciprocal distance to high achieving versus low achieving students (here gender would be ROWS and academic achievement would be COLS).
#' @param Option Whether to obtain the average ("Mean") or sum ("Sum") for the reciprocal distance matrix.
#'
#' @return A single value or matrix calculating the sum or average of compactness across all nodes, or split by category (or node).
#' @importFrom sna geodist
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' Carter T. Butts (2020). sna: Tools for Social Network Analysis. R package version 2.6. http://CRAN.R-project.org/package=sna
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xConnectedness()], [sna::geodist()]
#'
#' @examples
#' ## Two simple undirected one-mode networks (with the same value for connectedness)
#' xCompactness(ASNR_Fig10x3C)
#' xCompactness(ASNR_Fig10x3D)
#' ## Adding an attribute for rows and columns:
#' xCompactness(ASNR_Fig10x3D,ROWS=c(1,1,1,2,2,2),COLS=c(1,1,1,2,2,2))
#'
#' ## Some more complex examples with a two-mode network
#' NET2<-matrix(c(0,1,0,0,0,0,0, 1,0,0,0,0,0,0, 1,1,1,0,0,1,0, 0,1,1,0,0,0,0,
#'                0,0,0,1,1,0,0, 1,0,0,0,0,1,0, 0,1,1,1,0,1,0, 0,0,1,0,0,0,0),8,7,byrow=TRUE)
#' ## When running the function on the twomode, note the use of "A" and "B" to separate
#' ## mode A (row) from mode B (column) nodes:
#' xCompactness(NET2)
#' ## When adding an attribute for ROWS, note that mode A nodes are now split further
#' ## into 2 groups: A_1, A_2 and A_3.
#' xCompactness(NET2,ROWS=c(1,1,1,2,2,2,3,3))
#' ## We can also keep the values separate for each unique node in the column:
#' xCompactness(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS="Unique")
#' ## Or have a different attribute for mode A and mode B:
#' xCompactness(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS=c("P","E","T","T","P","T","E"))
#' ## Or combine unique (in the rows) with aggregating all (in the columns):
#' xCompactness(NET2,ROWS="Unique")
#'
xCompactness<-function(NET1, ROWS="All", COLS="All", Option="Mean")
{
  options(scipen = 100)
  if(!Option %in% c("Mean", "Sum")) {stop(' .  ### The argument (Option) should be "Mean" or "Sum" ###')}
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  if(sum(is.na(NET1))>0) {stop(' .  ## Matrix should not contain NAs')}
  if(sum(NET1!=0 & NET1!=1)>0) {stop(' .  ## Matrix should be binary (0 and 1 values only)')}

  #NETWORK MATRIX CHECK
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  #CALCULATION
  GDISTNI<-1/((sna::geodist(NET1))$gdist)

  #SPECIAL WAY TO DEAL WITH TWO-MODE -> BIPARTITE
  if(NR1!=NC1)
  {
    if(length(ROWS)==1)
    {
      if(ROWS=="All") {AG_R<-1}
      if(ROWS=="Unique") {AG_R<-(-1)}
    }
    else {AG_R<-0}

    if(length(COLS)==1)
    {
      if(COLS=="All") {AG_C<-1}
      if(COLS=="Unique") {AG_C<-(-1)}
    }
    else {AG_C<-0}
    # NOW ADD A AND BE BEFORE
    if(AG_R==1){PART_A<-rep("A",NR1)}
    if(AG_R==-1 & is.null(rownames(NET1))) {PART_A<-paste("A",c(1:NR1),sep="_a")}
    if(AG_R==-1 & !is.null(rownames(NET1))) {PART_A<-paste("A",rownames(NET1),sep="_")}
    if(AG_R==0) {PART_A<-paste("A",ROWS,sep="_")}

    if(AG_C==1){PART_B<-rep("B",NC1)}
    if(AG_C==-1 & is.null(colnames(NET1))) {PART_B<-paste("B",c(1:NC1),sep="_e")}
    if(AG_C==-1 & !is.null(colnames(NET1))) {PART_B<-paste("B",colnames(NET1),sep="_")}
    if(AG_C==0) {PART_B<-paste("B",COLS,sep="_")}

    ROWS<-c(PART_A,PART_B)
    COLS<-c(PART_A,PART_B)
  }
  OUTPUT1<-xDensity(GDISTNI, Loops=FALSE, ROWS=ROWS, COLS=COLS, Option=Option)
  OUTPUT1
}


#--------------------------------------------------------------------
#' Calculates the proportion (or number) of dyads that can reach each other in k steps in a network or in a subparts of a network.
#'
#' Using a one-mode or two-mode binary network, this function calculates the proportion (or number) of dyads that reach each other in a maximum of k steps.
#' Nodes that can reach each other in k steps are given a value 1, while those that cannot reach each other in k steps get a value 0 in the kreach matrix.
#' In practice it takes the geodesic distance matrix, and dichotomized this matrix using kreach, and then takes the average (Option="Mean") or the sum (Option="Sum") of the values for all pairs in the network.
#' - For a two-mode network, the result can be split into three parts (within mode A, within mode B, and between mode A and mode B nodes), and therefore we will obtain a matrix with within and between values.
#' - Missing data NA are not allowed.
#' - In addition, if an attribute file is specified for rows (ROWS) and/or the columns (COLS), then the average or sum is provided for each unique group in that attribute file. This might be useful to examine whether specific types of actors tend on average to be able to reach specific other nodes in k steps. Note that ROWS and COLS does not have to contain the same attribute.
#' - If ROWS or COLS is set to "Unique" (default is "All"), then the matrix is not aggregated for that row or column, and so the k-reach to other (sets of) nodes is obtained for each node. This might be useful if you want to know, for example, if specific nodes tend to be within k steps from others nodes with specific types.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#' @param kReach A single cut-off value.
#' @param ROWS Whether to calculate the sum or average of the k-reach over all nodes in the row (="All"), whether to keep each node in the row as a unique entity ("Unique"), or whether to aggregate based on some attribute variable (in which case (ROWS= name of the attribute dataset).
#' @param COLS Same as for rows, but here for columns. Note that ROWS and COLS can contain a different attribute variable, for example to see whether there is any difference in male or female students in their tendency to reach high achieving versus low achieving students within k steps (here gender would be ROWS and academic achievement would be COLS).
#' @param Option Whether to obtain the average ("Mean") or sum ("Sum") for the reciprocal distance matrix.
#'
#' @return A single value or matrix calculating the sum or average of compactness across all nodes, or split by category (or node).
#' @importFrom sna geodist
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' Carter T. Butts (2020). sna: Tools for Social Network Analysis. R package version 2.6. http://CRAN.R-project.org/package=sna
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xConnectedness()], [xUCINET::xCompactness()], [sna::geodist()]
#'
#' @examples
#' ## Two simple undirected one-mode networks (with the same value for connectedness)
#' xProportionReach(ASNR_Fig10x3C, kReach=2)
#' xProportionReach(ASNR_Fig10x3D, kReach=3)
#' ## Adding an attribute for rows and columns:
#' xProportionReach(ASNR_Fig10x3D,ROWS=c(1,1,1,2,2,2),COLS=c(1,1,1,2,2,2))
#'
#' ## Some more complex examples with a two-mode network
#' NET2<-matrix(c(0,1,0,0,0,0,0, 1,0,0,0,0,0,0, 1,1,1,0,0,1,0, 0,1,1,0,0,0,0,
#'                0,0,0,1,1,0,0, 1,0,0,0,0,1,0, 0,1,1,1,0,1,0, 0,0,1,0,0,0,0),8,7,byrow=TRUE)
#' ## When running the function on the twomode, note the use of "A" and "B" to separate
#' ## mode A (row) from mode B (column) nodes:
#' xProportionReach(NET2)
#' ## When adding an attribute for ROWS, note that mode A nodes are now split further
#' ## into 2 groups: A_1, A_2 and A_3.
#' xProportionReach(NET2,ROWS=c(1,1,1,2,2,2,3,3))
#' ## We can also keep the values separate for each unique node in the column:
#' xProportionReach(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS="Unique")
#' ## Or have a different attribute for mode A and mode B:
#' xProportionReach(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS=c("P","E","T","T","P","T","E"))
#' ## Or combine unique (in the rows) with aggregating all (in the columns):
#' xProportionReach(NET2,ROWS="Unique")
#'
xProportionReach<-function(NET1, kReach=2, ROWS="All", COLS="All", Option="Mean")
{
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  if(!Option %in% c("Mean", "Sum")) {stop(' .  ### The argument (Option) should be "Mean" or "Sum" ###')}
  if(sum(is.na(NET1))>0) {stop(' .  ## Matrix should not contain NAs')}
  if(sum(NET1!=0 & NET1!=1)>0) {stop(' .  ## Matrix should be binary (0 and 1 values only)')}
  #NETWORK MATRIX CHECK
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  #CALCULATION
  GDISTNI<-((sna::geodist(NET1))$gdist)<=kReach

  #SPECIAL WAY TO DEAL WITH TWO-MODE -> BIPARTITE
  if(NR1!=NC1)
  {
    if(length(ROWS)==1)
    {
      if(ROWS=="All") {AG_R<-1}
      if(ROWS=="Unique") {AG_R<-(-1)}
    }
    else {AG_R<-0}

    if(length(COLS)==1)
    {
      if(COLS=="All") {AG_C<-1}
      if(COLS=="Unique") {AG_C<-(-1)}
    }
    else {AG_C<-0}
    # NOW ADD A AND BE BEFORE
    if(AG_R==1){PART_A<-rep("A",NR1)}
    if(AG_R==-1 & is.null(rownames(NET1))) {PART_A<-paste("A",c(1:NR1),sep="_a")}
    if(AG_R==-1 & !is.null(rownames(NET1))) {PART_A<-paste("A",rownames(NET1),sep="_")}
    if(AG_R==0) {PART_A<-paste("A",ROWS,sep="_")}

    if(AG_C==1){PART_B<-rep("B",NC1)}
    if(AG_C==-1 & is.null(colnames(NET1))) {PART_B<-paste("B",c(1:NC1),sep="_e")}
    if(AG_C==-1 & !is.null(colnames(NET1))) {PART_B<-paste("B",colnames(NET1),sep="_")}
    if(AG_C==0) {PART_B<-paste("B",COLS,sep="_")}

    ROWS<-c(PART_A,PART_B)
    COLS<-c(PART_A,PART_B)
  }
  OUTPUT1<-xDensity(GDISTNI, Loops=FALSE, ROWS=ROWS, COLS=COLS, Option=Option)
  OUTPUT1
}


#--------------------------------------------------------------------
#' Calculates the (out)degree centralization for a one-mode or two-mode network.
#'
#' Using a one-mode or two-mode binary network, this function calculates the level of centralization in a network based on (out)degree.
#' It takes the maximum observed (out)degree as a point of reference against which the (out)degrees for all other nodes are compared.
#' No missing data allowed. The Loops is ignored.
#'
#' @param NET1 A one-mode or two-mode network stored as a 'matrix' object.
#'
#' @return A single value.
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xDegreeCentrality()]
#'
#' @examples
#' ## Two simple undirected one-mode networks (with the same value for connectedness)
#' xProportionReach(ASNR_Fig10x3C, kReach=2)
#' xProportionReach(ASNR_Fig10x3D, kReach=3)
#' ## Adding an attribute for rows and columns:
#' xProportionReach(ASNR_Fig10x3D,ROWS=c(1,1,1,2,2,2),COLS=c(1,1,1,2,2,2))
#'
#' ## Some more complex examples with a two-mode network
#' NET2<-matrix(c(0,1,0,0,0,0,0, 1,0,0,0,0,0,0, 1,1,1,0,0,1,0, 0,1,1,0,0,0,0,
#'                0,0,0,1,1,0,0, 1,0,0,0,0,1,0, 0,1,1,1,0,1,0, 0,0,1,0,0,0,0),8,7,byrow=TRUE)
#' ## When running the function on the twomode, note the use of "A" and "B" to separate
#' ## mode A (row) from mode B (column) nodes:
#' xProportionReach(NET2)
#' ## When adding an attribute for ROWS, note that mode A nodes are now split further
#' ## into 2 groups: A_1, A_2 and A_3.
#' xProportionReach(NET2,ROWS=c(1,1,1,2,2,2,3,3))
#' ## We can also keep the values separate for each unique node in the column:
#' xProportionReach(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS="Unique")
#' ## Or have a different attribute for mode A and mode B:
#' xProportionReach(NET2,ROWS=c(1,1,1,2,2,2,3,3),COLS=c("P","E","T","T","P","T","E"))
#' ## Or combine unique (in the rows) with aggregating all (in the columns):
#' xProportionReach(NET2,ROWS="Unique")
#'
xCentralization<-function(NET1)
{
  options(scipen = 100)
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  if(sum(is.na(NET1))>0) {stop(' .  ## Matrix should not contain NAs')}
  if(sum(NET1!=0 & NET1!=1)>0) {stop(' .  ## Matrix should be binary (0 and 1 values only)')}
  #NETWORK MATRIX CHECK
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(NR1==NC1){if(NET1==t(NET1)){NV<-2}else{NV<-1}}else{NV<-0}
  if(NR1==NC1){diag(NET1)<-0}
  #CALCULATION
  DCEN1<-rowSums(NET1,na.rm=T)
  OUTPUT1<-sum(max(DCEN1)-DCEN1)/((NR1-1)*(NC1-NV)) # if directed (NR1-1)*(NC-1), if undirected needs to be (NR1-1)*(NC1-2), if two mode (NR1-1)*(NC)
  OUTPUT1
}

#--------------------------------------------------------------------
#' Calculates measures of homophily for a one-mode network.
#'
#' Using a one-mode binary network, this function calculates measures of homophily.
#' This includes the proportion homophilous ties, the EI index and Yule's Q.
#' Missing data are allowed. The Loops is ignored.
#' Provides output:
#'  -"a": the number of pairs of nodes that have a tie and where both nodes belong to the same category.
#'  -"b": the number of pairs of nodes that have a tie and where both nodes belong to a different category.
#'  -"c": the number of pairs of nodes that do not have a tie and where both nodes belong to the same category.
#'  -"d": the number of pairs of nodes that do not have a tie and where both nodes belong to a different category.
#   -"PctSame": Percentage of all pairs of connected nodes that belong to the same type as ego (a/(a+b)).
#   -"EIIndex": Classic measure that takes the number of pairs of connected nodes that belong to a different
#      type (External = b), and subtracts the number of pairs of connected nodes that belong to the same type (Internal = a).
#      This difference is divided by the total number of connected pairs (a+b).
#   -"YulesQ": [(a*d)-(b*c)]/[(a*d)+(b*c)].
#'
#' @param NET1 A one-mode binary network stored as a 'matrix' object.
#' @param ATTR1 A vector containing the attribute.
#'
#' @return A vector with measures related to homophily.
#' @export
#'
#' @references
#' Chapter 10. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#'
#' @seealso [xUCINET::xDensity()], [xUCINET::xEgoAlterSimilarityCat()]
#'
#' @examples
#' ## Two examples of directed networks (See chapter 10 of Borgatti et al., 2022):
#' xGroupHomophily(Krackhardt_HighTech$Advice,Krackhardt_HighTech$Attributes$Department)
#' xGroupHomophily(Krackhardt_HighTech$Friendship,Krackhardt_HighTech$Attributes$Department)
#'
xGroupHomophily<-function(NET1,ATTR1)
{
  options(scipen = 100)
  #NETWORK MATRIX CHECK
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(NR1!=NC1) {stop(' .  ### Network file NET1 needs to be a one-mode network ###')}
  diag(NET1)<-NA

  #---------------------------------------------------------------------------------------------------
  #ATTRIBUTES CHECK
  if(!is.matrix(ATTR1) & !is.vector(ATTR1)) {stop('\n',' .  ###  Attribute file (ATTR1) needs to be of class "matrix" or "vector"','\n')}
  if(is.matrix(ATTR1)){
    if(min(dim(ATTR1))!=1) {stop('\n',' .  ###  Attribute file (ATTR1) needs to be a single row or column','\n')}
  }
  if(is.matrix(ATTR1)) {NAR<-max(dim(ATTR1))} else {NAR<-length(ATTR1)}
  if(NAR!=NR1) {stop('\n',' .  ###  Attribute file (ATTR1) needs to have the same length as the number of rows/columns of the matrix (NET1)','\n')}

  #CONSIDER ALL POSSIBILITIES:
  SUM1<-sapply(by(NET1,ATTR1,FUN=colSums, na.rm=T),identity)
  SUMF<-sapply(by(SUM1,ATTR1,FUN=colSums, na.rm=T),identity)

  SUM1VAL<-sapply(by(!is.na(NET1),ATTR1,FUN=colSums, na.rm=T),identity)
  SUMFVAL<-sapply(by(SUM1VAL,ATTR1,FUN=colSums, na.rm=T),identity)
  #now create a,b,c,d
  OUTPUT1<-c(sum(diag(SUMF)),
             sum(SUMF)-sum(diag(SUMF)),
             sum(diag(SUMFVAL))-sum(diag(SUMF)),
             sum(SUMFVAL)-sum(diag(SUMFVAL))-(sum(SUMF)-sum(diag(SUMF))))
  #add indices
  OUTPUT1<-matrix(c(
    OUTPUT1,
    OUTPUT1[1]/(OUTPUT1[1]+OUTPUT1[2]),
    (OUTPUT1[2]-OUTPUT1[1])/(OUTPUT1[2]+OUTPUT1[1]),
    (OUTPUT1[4]+OUTPUT1[2]-OUTPUT1[3]-OUTPUT1[1])/(sum(!is.na(NET1))),
    ((OUTPUT1[1]*OUTPUT1[4])-(OUTPUT1[2]*OUTPUT1[3]))/((OUTPUT1[1]*OUTPUT1[4])+(OUTPUT1[2]*OUTPUT1[3]))
  ),
  8,1)
  rownames(OUTPUT1)<-c("a","b","c","d","%homoph","EI","ExpectEI","YulesQ")
  OUTPUT1
}
