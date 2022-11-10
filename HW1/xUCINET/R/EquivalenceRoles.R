#' Calculate structural equivalence for one or more networks
#'
#' Computes for each pair of nodes (i and j) the level of structural equivalence based upon comparisons of rows
#' (and potentially also columns) in the adjacency matrix (or whatever else is used as input matrix).
#' Structural equivalence is based on a pair of nodes (i and j) being in the same/in a similar way connected to
#' other nodes (k).
#' Multiple relations are permissible as input (which should be combined in a list) and in such a case the
#' combined similarity/distance across these different relations are considered. The matrix can contain values or
#' binary data.
#' When consider a pair of nodes (i and j), the measure can include or exclude reciprocal connections (i->j and j->i)
#' as well as selfloops (i->i) and (j->j).
#' Calculations of similarity can be based on a a variety of measures, including Pearson correlation, exact matches
#' or matches of positive matches only.
#' Note that measures based on absolute difference and Euclidean distances produces a distance matrix, while all
#' the other options produce a similarity matrix.
#'
#' @param NET1 A matrix or list of matrices, which could be part of a xUCINET project
#' @param Method Method used to obtain the similarity/distances between two nodes, i and j, which can be based on:
#'        \itemize{
#'         \item "AbsDiff" the sum of the absolute differences between both nodes in their relations to others. I.e.,
#'         the value for the pair of nodes i,j is given by: |(i->k)-(j->k)| taken over all k. In the case where also
#'         the incoming ties are considered (see IncludeTransposed) this also includes: |(k->i)-(k->j)| taken over
#'         all k. The output is a distance matrix.
#'         \item "Euclidean" the squared square root of the sum of the squares of the differences between
#'         corresponding values. I.e., the sum of the value for the pair of nodes i,j is given by: ((i->k)-(j->k))^2
#'         taken over all k (as well as ((k->i)-(k->j))^2 in the case when also incoming ties are considered).
#'         For the resulting sum for each pair i,j the square root is then taken. The output is a distance matrix.
#'         \item "MatchesN" the number of exact matches, i.e., the number of times that i and j relate to other nodes
#'         k with exactly the same value. (i->k)==(j->k) taken over all k (as well as (k->i)==(k->j) if the incoming
#'         ties are also considered). The output is a similarity matrix.
#'         \item "PosMatchesN" the number of exact positive matches, i.e., the number of times that i and j relate to
#'         other nodes k with exactly the same positive value (i.e., value different from 0). This approach differs
#'         from "MatchesN" in that it excludes similarity based on two nodes not being connected to others. The output
#'         is a similarity matrix. For a nimary network it provides an index of the number of common non-zero choices.
#'         \item  "Product" the product of the link between i and k, and j and k. I.e., (i->k)*(j->k) taken over all k
#'         (and (k->i)-(k->j) if the incoming ties are also considered). The output is a similarity matrix with high
#'         values indicating both i and j are highly connected to k.
#'         \item  "Pearson" takes the values for the links between i and k, and correlates it with the values for the
#'         corresponding links between j and k. I.e., cor ((i->k),(j->k)) for outgoing ties.  The output is a similarity
#'          matrix.
#'         \item  "Kendall" similar to "Pearson", except it uses Kendall's Tau.
#'         \item  "Spearman" similar to "Pearson", except it uses Spearman's Rho, i.e., the ranked value.
#'         }
#' @param IncludeTransposed Whether or not to focus only on outgoing ties to define the structural equivalence
#'         between i and j (based on the measure defined in Method. If FALSE then only the outgoing ties to other
#'         nodes (k) are considered (i.e., i->k and j->k) when calculating the structural equivalence between nodes
#'         i and j, while if TRUE the also the incoming ties to i and j from other nodes (i.e., k->i and k->j) will
#'         be used to calculate the structural equivalence between nodes i and j.
#' @param Choiceij What to do with ties between i and j and selfloops?
#'        \itemize{
#'         \item If "Ignore", then only the ties from i and j with other nodes k are considered (where k is
#'         not i or j), i.e., the values for i->i, j->i, i->j and j->j are ignored.
#'         \item If "OnlyReciprocal", then loops (i->i and j->j) are ignored, and i->j is compared to j->i.
#'         \item If "Reciprocal", then i->j is compared to j->i and i->i with j->j.
#'         \item If "Original", then the classic approach of comparing rows/columns is used, and hence i->i is
#'         compared to j->i, while i->j is compared to j->j (plus i->i with i->j and j->i with j->j when incoming
#'         ties are also considering)
#'         }

#' @param DigitRound Whether to round the output to a specified number of decimal places (default = 5).
#'
#' @return A matrix containing the level of similarity between each pair of nodes (or distances in the case of
#' choosing the euclidean distance or absolute difference).
#'
#' @importFrom utils combn
#' @importFrom stats cor
#'
#' @export
#'
#' @references
#' Chapter 12. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @examples
#' ## A small example with 3 networks focusing on outgoing ties only (default):
#' xStructuralEquivalence(ASNR_Fig12x1R1)
#' xStructuralEquivalence(list(ASNR_Fig12x1R1,ASNR_Fig12x1R2,ASNR_Fig12x1R3))
#'
#' ## Incoming and outgoing ties:
#' xStructuralEquivalence(list(ASNR_Fig12x1R1,ASNR_Fig12x1R2,ASNR_Fig12x1R3))
#'
#' ## Apply the function using both incoming and outgoing ties for esteem and disesteem
#' ## (Sampson Monastery):
#' Samp1<-list(Sampson_Monastery$Esteem, Sampson_Monastery$Disesteem)
#' xStructuralEquivalence(Samp1,IncludeTransposed=TRUE)

xStructuralEquivalence<-function(NET1, Method="Euclidean",
                                 IncludeTransposed=F, Choiceij="OnlyReciprocal",
                                 DigitRound=5)
{
  # -------------------------------------------------------------
  #check input and change list into array so we can easily extract parts [sender,receiver,matrixnumber]
  if (is.list(NET1)==F & is.matrix(NET1)==F)
  {
    stop('NET1 object needs to be either of class "matrix" or a list of matrices')
  }

  if (is.list(NET1)==T)
  {
    if(sum(dim(NET1[[1]])!=unlist(lapply(NET1, dim)))>0)
    {
      stop('NET1 object needs to be a list where each element is a matrix with the same dimensions')
    }
  }

  TwoMode<-FALSE

  if (is.matrix(NET1)==T)
  {
    if (dim(NET1)[1]!=dim(NET1)[2])
    {
      TwoMode<-TRUE
      IncludeTransposed<-FALSE
      Choiceij="Original"
    }
    ARRAY1<-array(NET1,dim=c(dim(NET1),1))
    ROWNAMES<-rownames(NET1)
  }

  if (is.list(NET1)==T)
  {
    if (dim(NET1[[1]])[1]!=dim(NET1[[1]])[2])
    {
      TwoMode<-TRUE
      IncludeTransposed<-FALSE
      Choiceij="Original"
    }
    ARRAY1<-array(unlist(NET1), dim=c(dim(NET1[[1]]),length(NET1)) )
    ROWNAMES<-rownames(NET1[[1]])
  }

  #==== 1.1. Check all the measures to be included from the full list below ====
  MethodSEAll<-c("AbsDiff","Euclidean","MatchesN","PosMatchesN","Product","Pearson","Kendall","Spearman")
  #Check if the Method is a valid option
  if ((sum(Method %in% MethodSEAll))!=1)
  {
    stop("WARNING: The 'Method' could not be identified, and might contain a typo:", "\n")
  }

  ChoiceijAll<-c("Ignore", "OnlyReciprocal", "Reciprocal", "Original")
  #Check if the Choiceij is a valid option
  if ((sum(Choiceij %in% ChoiceijAll))!=1)
  {
    stop("WARNING: The 'Choiceij' could not be identified, and might contain a typo:", "\n")
  }

  if(TwoMode==FALSE)
  {
    #-----------------------
    # Data
    # Change diagonal into NA, so any analysis with i->i, and j->j are ignored when using matrices.
    # First copy the data to a matrix called ARRAY1 and another called ARRAY2:
    ARRAY2<-ARRAY1
    Nmat<-dim(ARRAY1)[3]
    #Put diagonal to NA for matrices
    for(k in 1:Nmat)
    {
      diag(ARRAY1[,,k])<-NA
    }

    # -------------------------------------------------------------
    # Decision which Method is needed
    Nact<-dim(ARRAY1)[1]

    # CHOICE 1. Euclidean distance ====
    if(Method=="Euclidean")
    {

      # Set up the OUTPUT1 matrix:
      OUTPUT1<-matrix(0,Nact,Nact)

      # CHOICE 1.1.
      # If Choiceij=="Ignore", then all 4 cells are ignored (i->i, j->i, i->j and j->j)
      if(Choiceij=="Ignore")
      {

        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,])-c(ARRAY1[x[2],,]))^2 )),na.rm=T)^.5
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[,x[1],])-c(ARRAY1[x[2],,],ARRAY1[,x[2],]))^2 )),na.rm=T)^.5
        }
      }

      # CHOICE 1.2.
      # If "OnlyReciprocal", then reciprocal i->j with j->i (plus j->i with i->j when also considering incoming) with i->i and j->j ignored
      if(Choiceij=="OnlyReciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],])-c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],]))^2 )),na.rm=T)^.5
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],])-c(ARRAY1[x[2],,],ARRAY1[,x[2],],ARRAY1[x[2],x[1],],ARRAY1[x[1],x[2],]))^2 )),na.rm=T)^.5
        }
      }

      # CHOICE 1.3.
      # If "Reciprocal", then reciprocal i->j with j->i and i->i with j->j (plus j->i with i->j and j->j with i->i when also considering incoming)
      if(Choiceij=="Reciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],])
                 -c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],]))^2 )),na.rm=T)^.5
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],])
                 -c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY1[x[1],x[2],],ARRAY2[x[2],x[2],]))^2 )),na.rm=T)^.5
        }
      }

      # CHOICE 1.4.
      # If "Original", then the classic approach of comparing rows/columns is used and compare i->i with j->i and i->j with j->j (plus i->i with i->j and j->i with j->j when also considering incoming)
      if(Choiceij=="Original")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],])
                 -c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],]))^2 )),na.rm=T)^.5
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],])
                 -c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY2[x[1],x[2],],ARRAY2[x[2],x[2],]))^2 )),na.rm=T)^.5
        }
      }

    }

    # CHOICE 2. Absolute Difference ====
    if(Method=="AbsDiff")
    {

      # Set up the OUTPUT1 matrix:
      OUTPUT1<-matrix(0,Nact,Nact)

      # CHOICE 2.1.
      # If Choiceij=="Ignore", then all 4 cells are ignored (i->i, j->i, i->j and j->j)
      if(Choiceij=="Ignore")
      {

        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,])-c(ARRAY1[x[2],,])) )   ),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY1[,x[1],])-c(ARRAY1[x[2],,],ARRAY1[,x[2],])))),na.rm=T)
        }
      }

      # CHOICE 2.2.
      # If "OnlyReciprocal", then reciprocal i->j with j->i (plus j->i with i->j when also considering incoming) with i->i and j->j ignored
      if(Choiceij=="OnlyReciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],])-c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],])-c(ARRAY1[x[2],,],ARRAY1[,x[2],],ARRAY1[x[2],x[1],],ARRAY1[x[1],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 2.3.
      # If "Reciprocal", then reciprocal i->j with j->i and i->i with j->j (plus j->i with i->j and j->j with i->i when also considering incoming)
      if(Choiceij=="Reciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],])
                    -c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],])
                    -c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY1[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 2.4.
      # If "Original", then the classic approach of comparing rows/columns is used and compare i->i with j->i and i->j with j->j (plus i->i with i->j and j->i with j->j when also considering incoming)
      if(Choiceij=="Original")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],])
                    -c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              (abs (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],])
                    -c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY2[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

    }

    # CHOICE 3. Exact Matches ====
    if(Method=="MatchesN")
    {

      # Set up the OUTPUT1 matrix:
      OUTPUT1<-matrix(0,Nact,Nact)

      # CHOICE 3.1.
      # If Choiceij=="Ignore", then all 4 cells are ignored (i->i, j->i, i->j and j->j)
      if(Choiceij=="Ignore")
      {

        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,]))==(c(ARRAY1[x[2],,])) )   ),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[,x[1],]))==(c(ARRAY1[x[2],,],ARRAY1[,x[2],])))),na.rm=T)
        }
      }

      # CHOICE 3.2.
      # If "OnlyReciprocal", then reciprocal i->j with j->i (plus j->i with i->j when also considering incoming) with i->i and j->j ignored
      if(Choiceij=="OnlyReciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],]))==(c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],]))==(c(ARRAY1[x[2],,],ARRAY1[,x[2],],ARRAY1[x[2],x[1],],ARRAY1[x[1],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 3.3.
      # If "Reciprocal", then reciprocal i->j with j->i and i->i with j->j (plus j->i with i->j and j->j with i->i when also considering incoming)
      if(Choiceij=="Reciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],]))
                ==(c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],]))
                ==(c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY1[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 3.4.
      # If "Original", then the classic approach of comparing rows/columns is used and compare i->i with j->i and i->j with j->j (plus i->i with i->j and j->i with j->j when also considering incoming)
      if(Choiceij=="Original")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],]))
                ==(c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],]))
                ==(c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY2[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

    }

    # CHOICE 4. Positive Matches ====
    if(Method=="PosMatchesN")
    {

      # Set up the OUTPUT1 matrix:
      OUTPUT1<-matrix(0,Nact,Nact)

      # CHOICE $.1.
      # If Choiceij=="Ignore", then all 4 cells are ignored (i->i, j->i, i->j and j->j)
      if(Choiceij=="Ignore")
      {

        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,])>0) * (c(ARRAY1[x[1],,])==c(ARRAY1[x[2],,])) )   ),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY1[,x[1],])>0) * (c(ARRAY1[x[1],,],ARRAY1[,x[1],])==c(ARRAY1[x[2],,],ARRAY1[,x[2],])))),na.rm=T)
        }
      }

      # CHOICE 4.2.
      # If "OnlyReciprocal", then reciprocal i->j with j->i (plus j->i with i->j when also considering incoming) with i->i and j->j ignored
      if(Choiceij=="OnlyReciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],])>0)* (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],])==c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],])>0)* (c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],])==c(ARRAY1[x[2],,],ARRAY1[,x[2],],ARRAY1[x[2],x[1],],ARRAY1[x[1],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 4.3.
      # If "Reciprocal", then reciprocal i->j with j->i and i->i with j->j (plus j->i with i->j and j->j with i->i when also considering incoming)
      if(Choiceij=="Reciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],])>0)* (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],])
                                                                             ==c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],])>0)* (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],])
                                                                                                                                  ==c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY1[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 4.4.
      # If "Original", then the classic approach of comparing rows/columns is used and compare i->i with j->i and i->j with j->j (plus i->i with i->j and j->i with j->j when also considering incoming)
      if(Choiceij=="Original")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],])>0) * (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],])
                                                                              ==c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ((c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],])>0) * (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],])
                                                                                                                                   ==c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY2[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

    }

    # CHOICE 5. Product ====
    if(Method=="Product")
    {
      # Set up the OUTPUT1 matrix:
      OUTPUT1<-matrix(0,Nact,Nact)

      # CHOICE 5.1.
      # If Choiceij=="Ignore", then all 4 cells are ignored (i->i, j->i, i->j and j->j)
      if(Choiceij=="Ignore")
      {

        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,]))*(c(ARRAY1[x[2],,])) )   ),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[,x[1],]))*(c(ARRAY1[x[2],,],ARRAY1[,x[2],])))),na.rm=T)
        }
      }

      # CHOICE 5.2.
      # If "OnlyReciprocal", then reciprocal i->j with j->i (plus j->i with i->j when also considering incoming) with i->i and j->j ignored
      if(Choiceij=="OnlyReciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],]))*(c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],]))*(c(ARRAY1[x[2],,],ARRAY1[,x[2],],ARRAY1[x[2],x[1],],ARRAY1[x[1],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 5.3.
      # If "Reciprocal", then reciprocal i->j with j->i and i->i with j->j (plus j->i with i->j and j->j with i->i when also considering incoming)
      if(Choiceij=="Reciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],]))
                *(c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],]))
                *(c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY1[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

      # CHOICE 5.4.
      # If "Original", then the classic approach of comparing rows/columns is used and compare i->i with j->i and i->j with j->j (plus i->i with i->j and j->i with j->j when also considering incoming)
      if(Choiceij=="Original")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],]))
                *(c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            colSums(combn(Nact,2,FUN=function(x)
              ( (c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],]))
                *(c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY2[x[1],x[2],],ARRAY2[x[2],x[2],])) )),na.rm=T)
        }
      }

    }

    # CHOICE 6. Correlation ====
    if(Method=="Pearson"|Method=="Kendall"|Method=="Spearman")
    {
      # Set up the OUTPUT1 matrix:
      OUTPUT1<-matrix(1,Nact,Nact)

      #Make capital letters small to fit cor function
      Method1<-tolower(Method)

      # CHOICE 6.1.
      # If Choiceij=="Ignore", then all 4 cells are ignored (i->i, j->i, i->j and j->j)
      if(Choiceij=="Ignore")
      {

        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (stats::cor(
                c(ARRAY1[x[1],,]),c(ARRAY1[x[2],,])
                ,use="pairwise.complete.obs", method=Method1) ))
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY1[,x[1],]),c(ARRAY1[x[2],,],ARRAY1[,x[2],])
                ,use="pairwise.complete.obs", method=Method1) ))

        }
      }

      # CHOICE 6.2.
      # If "OnlyReciprocal", then reciprocal i->j with j->i (plus j->i with i->j when also considering incoming) with i->i and j->j ignored
      if(Choiceij=="OnlyReciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],]),c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],])
                ,use="pairwise.complete.obs", method=Method1) ))
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY1[,x[1],],ARRAY1[x[1],x[2],],ARRAY1[x[2],x[1],]),
                c(ARRAY1[x[2],,],ARRAY1[,x[2],],ARRAY1[x[2],x[1],],ARRAY1[x[1],x[2],])
                ,use="pairwise.complete.obs", method=Method1) ))
        }
      }

      # CHOICE 6.3.
      # If "Reciprocal", then reciprocal i->j with j->i and i->i with j->j (plus j->i with i->j and j->j with i->i when also considering incoming)
      if(Choiceij=="Reciprocal")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],]),
                c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],])
                ,use="pairwise.complete.obs", method=Method1) ))
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY1[x[1],x[2],],ARRAY2[x[1],x[1],],ARRAY1[,x[1],],ARRAY1[x[2],x[1],],ARRAY2[x[1],x[1],]),
                c(ARRAY1[x[2],,],ARRAY1[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY1[x[1],x[2],],ARRAY2[x[2],x[2],])
                ,use="pairwise.complete.obs", method=Method1) ))
        }
      }

      # CHOICE 6.4.
      # If "Original", then the classic approach of comparing rows/columns is used and compare i->i with j->i and i->j with j->j (plus i->i with i->j and j->i with j->j when also considering incoming)
      if(Choiceij=="Original")
      {
        # DECIDE ON OUTGOING OR BOTH OUTGOING AND INCOMING?
        if(IncludeTransposed==F)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],]),
                c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],])
                ,use="pairwise.complete.obs", method=Method1) ))
        }

        if(IncludeTransposed==T)
        {
          OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
            combn(Nact,2,FUN=function(x)
              (cor(
                c(ARRAY1[x[1],,],ARRAY2[x[1],x[1],],ARRAY2[x[1],x[2],],ARRAY1[,x[1],],ARRAY2[x[1],x[1],],ARRAY2[x[2],x[1],]),
                c(ARRAY1[x[2],,],ARRAY2[x[2],x[1],],ARRAY2[x[2],x[2],],ARRAY1[,x[2],],ARRAY2[x[1],x[2],],ARRAY2[x[2],x[2],])
                ,use="pairwise.complete.obs", method=Method1) ))
        }
      }

    }
  }


  if(TwoMode==TRUE)
  {
    #-----------------------
    ARRAY2<-ARRAY1
    Nmat<-dim(ARRAY1)[3]
    # -------------------------------------------------------------
    # Decision which Method is needed
    Nact<-dim(ARRAY1)[1]
    OUTPUT1<-matrix(0,Nact,Nact)

    # CHOICE 1. Euclidean distance ====
    if(Method=="Euclidean")
    {
      OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
        colSums(combn(Nact,2,FUN=function(x)
          ( (c(ARRAY1[x[1],,])-c(ARRAY1[x[2],,]))^2 )),na.rm=T)^.5
    }

    # CHOICE 2. Absolute Difference ====
    if(Method=="AbsDiff")
    {
      OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
        colSums(combn(Nact,2,FUN=function(x)
          (abs (c(ARRAY1[x[1],,])-c(ARRAY1[x[2],,])) )   ),na.rm=T)
    }

    # CHOICE 3. Exact Matches ====
    if(Method=="MatchesN")
    {
      OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
        colSums(combn(Nact,2,FUN=function(x)
          ( (c(ARRAY1[x[1],,]))==(c(ARRAY1[x[2],,])) )   ),na.rm=T)
    }

    # CHOICE 4. Positive Matches ====
    if(Method=="PosMatchesN")
    {
      OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
        colSums(combn(Nact,2,FUN=function(x)
          ((c(ARRAY1[x[1],,])>0) * (c(ARRAY1[x[1],,])==c(ARRAY1[x[2],,])) )   ),na.rm=T)
    }

    # CHOICE 5. Product ====
    if(Method=="Product")
    {
      OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
        colSums(combn(Nact,2,FUN=function(x)
          ( (c(ARRAY1[x[1],,]))*(c(ARRAY1[x[2],,])) )   ),na.rm=T)
    }

    # CHOICE 6. Correlation ====
    if(Method=="Pearson"|Method=="Kendall"|Method=="Spearman")
    {
      OUTPUT1[lower.tri(OUTPUT1,diag=F)]<-
        combn(Nact,2,FUN=function(x)
          (cor(
            c(ARRAY1[x[1],,]),c(ARRAY1[x[2],,])
            ,use="pairwise.complete.obs", method=Method1) ))
    }
  }


  #Now generate output (ensure to fill the other side of matrix and diagonal is correct)
  OUTPUT1<-OUTPUT1+t(OUTPUT1)
  rownames(OUTPUT1)<-ROWNAMES
  colnames(OUTPUT1)<-ROWNAMES
  cat('Note that "Euclidean" and "AbsDiff" result in a matrix where high values represent differences.')
  cat(' "MatchesN", "PosMatchesN", "Product", "Pearson", "Kendall" and "Spearman" result in a matrix where high values represent similarities./n')
  round(OUTPUT1, digits=DigitRound)
}

#' Construct a blockmodel
#'
#' Takes a network and an attribute file containing partitions.
#'
#' @param NET1 A network
#' @param ATTR1 An attribute file containing partition to be performed.
#'
#' @return plot
#' @import blockmodeling
#' @export
#' @references
#' Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' Doreian, P., Batagelj, V., & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge: Cambridge University Press.
#' Žiberna, A. (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' Žiberna, A. (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' Žiberna, A. (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002

#' @examples
#' NET1<-matrix(c(0,1,0,1,0,
#'                0,0,0,0,1,
#'                1,1,0,1,1,
#'                0,0,0,0,1,
#'                1,1,1,0,0),5,5)
#' xBlockmodel(NET1,ATTR1=c(1,2,1,2,1))

xBlockmodel<-function(NET1, ATTR1)
{
  if(is.matrix(NET1)){blockmodeling::plot.mat(NET1,ATTR1)}
  if(is.list(NET1))
  {
    for (mk in 1:length(NET1))
    {
      blockmodeling::plot.mat(NET1[[mk]],ATTR1)
    }
  }
}

#' Finds optimal blockmodel
#'
#' Takes a network and finds the optimal solution(s) given the options requested.
#'
#' @param NET1 A network stored as a matrix
#' @param NOG Number of groups/clusters to consider.
#' @param BlockTypes Types of blocks allowed. See blockmodeling::optRandomParC.
#' @param Options Whether Binary or Ratio.
#' @param NTRIALS Number of trials considered to find the optimal solution.
#'
#' @return Vector of different groupings based on optimal solutions.
#' @importFrom blockmodeling optRandomParC
#' @importFrom graphics barplot
#' @export
#' @references
#' Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#' Doreian, P., Batagelj, V., & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge: Cambridge University Press.
#' Žiberna, A. (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' Žiberna, A. (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' Žiberna, A. (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002
#'
#' @examples
#' NET1<-matrix(c(0,1,0,1,0,
#'                0,0,0,0,1,
#'                1,1,0,1,1,
#'                0,0,0,0,1,
#'                1,1,1,0,0),5,5)
#' xBlockmodelOptimizing(NET1,NOG=2,BlockTypes=c("nul","com"),Options="Binary")

xBlockmodelOptimizing<-function(NET1, NOG, BlockTypes, Options, NTRIALS=50)
{
  if(is.matrix(NET1)==F) {stop(' .  ### Network file NET1 needs to be of class "matrix" ###')}
  NR1<-dim(NET1)[1]
  NC1<-dim(NET1)[2]
  if(sum(is.na(NET1)>0)){stop(' .   ### Network file NET1 cannot contain missing data. ')}

  if(Options=="Binary"){Options1="bin"}
  if(Options=="Ratio"){Options1="val"}
  OUTPUT0<-blockmodeling::optRandomParC(M=NET1, k=NOG, blocks=BlockTypes, rep=NTRIALS, approaches="bin")
  graphics::barplot(table(OUTPUT0$err),main="Number of Errors")
  OUTPUT1<-matrix(NA,NR1,length(OUTPUT0$best))
  for (k in 1:length(OUTPUT0$best))
  {
    SOLk<-OUTPUT0$best[[k]]
    OUTPUT1[,k]<-SOLk$clu
  }
  rownames(OUTPUT1)<-rownames(NET1)
  colnames(OUTPUT1)<-paste("GROUP_",c(1:length(OUTPUT0$best)),"_E",OUTPUT0$best$best1$err,sep="")
  OUTPUT1
}

#' Performs a REGE solution
#'
#' Takes a network as matrix or list of networks and finds the REGE solution.
#'
#' @param NET1 A network stored as a matrix or list of matrices
#'
#' @return A matrix representing regular equivalence values.
#' @importFrom blockmodeling REGE.for
#' @export
#' @references
#' Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#' Doreian, P., Batagelj, V., & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge: Cambridge University Press.
#' Žiberna, A. (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' Žiberna, A. (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' Žiberna, A. (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002
#'
#' @examples
#' NET1<-matrix(c(0,1,0,1,0,
#'                0,0,0,0,1,
#'                1,1,0,1,1,
#'                0,0,0,0,1,
#'                1,1,1,0,0),5,5)
#' xREGE(NET1)


xREGE<-function(NET1)
{
  if(is.list(NET1))
  {
    NR1<-dim(NET1[[1]])[1]
    NC1<-dim(NET1[[2]])[1]
    NE1<-length(NET1)
    INPUT1<-array(as.numeric(unlist(NET1)), c(NR1,NC1,NE1))
  } else {INPUT1<-NET1}
  OUTPUT0<-blockmodeling::REGE.for(M=INPUT1,E=1)
  OUTPUT1<-OUTPUT0$E
  rownames(OUTPUT1)<-rownames(NET1[[1]])
  colnames(OUTPUT1)<-colnames(NET1[[1]])
  OUTPUT1*100
}

#' Finds a core-periphery structure
#'
#' Takes a network and finds the core-periphery structure using the correlation between the values and a vector with 0s and 1s (concentration).
#'
#' @param NET1 A network stored as a matrix
#'
#' @return The eigenvector of each node, the correlation and an optimal solution.
#'
#' @importFrom igraph eigen_centrality graph.adjacency
#' @importFrom stats cor
#'
#' @export
#' @references
#' Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022).
#' Analyzing Social Networks Using R. SAGE.
#'
#' @examples
#' NET1<-matrix(c(0,1,0,1,0,
#'                0,0,0,0,1,
#'                1,1,0,1,1,
#'                0,0,0,0,1,
#'                1,1,1,0,0),5,5)
#' xCorePeriphery(NET1)

xCorePeriphery<-function(NET1)
{
  options(scipen=100)
  NR1<-dim(NET1)[1]
  NET1i<-igraph::graph.adjacency(t(NET1))
  OUTPUT1<-matrix(NA,NR1,5)
  EV<-igraph::eigen_centrality(NET1i,directed=TRUE)$vector
  OUTPUT1[,1]<-EV
  ORN<-factor(OUTPUT1[,1])
  OUTPUT1[,2]<-as.numeric(ORN)
  EVO<-EV[order(EV,decreasing=TRUE)]
  plot(EVO)
  OUTPUT1[,4]<-NR1
  for(k in 2:(max(OUTPUT1[,2])))
  {
    OUTPUT1[OUTPUT1[,2]==k,3]<-cor(OUTPUT1[,1],(OUTPUT1[,2]>=k))
    OUTPUT1[OUTPUT1[,2]==k,4]<-NR1-sum(OUTPUT1[,2]<k)
  }
  OUTPUT1[,5]<-1*(OUTPUT1[,4]<=OUTPUT1[which(OUTPUT1[,3]==max(OUTPUT1[,3],na.rm=TRUE)),4])
  colnames(OUTPUT1)<-c("Value","Order","Correlation","NCorewithK","MaxSolution")
  rownames(OUTPUT1)<-rownames(NET1)
  OUTPUT1
}

