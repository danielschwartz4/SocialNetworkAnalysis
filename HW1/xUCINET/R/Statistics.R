#' Performs a correlation test using a classic (and a permutation-based) significance test between two vector variables
#'
#' This function takes two variables (vectors of equal length) and returns the correlation coefficient("CorrCoef") and the classic (two-tailed) significance test ("ClassicSign_2tailed").
#' If NPERM (number of permutations) is set to an integer value higher than 0 (default is 1000), a permutations-based significance is also provided,
#' which will provide information about the correlations coefficients for the permuted variables (including the histogram). A higher value than 1000 is generally advisable. Both the one-tailed ("AsSmall" or "AsLarge")
#' and the two-tailed ("AsExtreme") test are given. Note that cases with missing values are allowed, but that this can impact the degrees of freedom for the correlations obtained for the permutations.
#' For example, suppose the same 20% of cases are missing for both variables (i.e., unit non-response), then permuting these missing values will mean that a lot more cases will have a missing for (at least) one of both variables.
#' - Last updated: 26 April 2022.
#'
#' @param VEC1 A first vector of values.
#' @param VEC2 A second vector of values, of equal length as VEC1.
#' @param NPERM A positive integer value indicating the number of permutations to be used for the statistical test. If NPERM is set to NULL, NA or 0, no permutation-based significance test is performed.
#'
#' @return Returns an object containing the correlation coefficient ($CorrCoeff) and standard significance test (ClassicSign_2tailed), and in the case a permutation test is performed, in addition a significance based on the permutation ("AsSmall", "AsLarge", and "AsExtreme").
#' @importFrom stats cor cor.test runif
#' @importFrom graphics hist
#' @export
#' @author Filip Agneessens <filipagneessens2@gmail.com>
#'
#' @examples
#' RESULS1<-xCorrelation(c(2,4,2,3,3,1),c(4,2,3,2,1,5), NPERM=5000)
#' RESULS1
#' RESULS2<-xCorrelation(c(2,3,2,1,2,3,2,1),c(4,2,3,2,4,2,3,2), NPERM=5000)
#' RESULS2
#' xCorrelation(c(1,4,2,1),c(4,2,3,3), NPERM=0)
#' xCorrelation(c(1,4,2,1),c(4,2,3,3), NPERM=NA)
#' xCorrelation(c(1,4,2,1),c(4,2,3,3), NPERM=NULL)

xCorrelation<-function(VEC1,VEC2,NPERM=1000)
{
  options(scipen=999)
  if(is.null(NPERM)) {NPERM<-0} # changed to 0 because NULL is not a value
  if(is.na(NPERM)) {NPERM<-0} # changed to 0 because NA is not a value
  if(!is.vector(VEC1))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "VEC1" is not a vector.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if(!is.vector(VEC2))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "VEC2" is not a vector.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if((ceiling(NPERM)!=floor(NPERM)) | (NPERM<0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The number of permutations ("NPERM") requested is not a positive integer value.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  NCASES1<-length(VEC1)
  NCASES2<-length(VEC2)

  if(NCASES1!=NCASES2)
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The length for VEC1 is not the same as the length for VEC2.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  COR1<-cor.test(VEC1,VEC2)
  CORC<-COR1$estimate
  names(CORC)<-NULL

  NVALID1<-sum(!is.na(VEC1))
  NVALID2<-sum(!is.na(VEC2))
  NVALID12<-sum(!is.na(VEC1*VEC2))

  # Only perform a permutation test if NPERM>0 (and not NULL or NA).
  if(NPERM>0)
  {
    OUTPUT0<-replicate(NPERM,cor(sample(VEC1),VEC2, use = "pairwise.complete.obs"))
    OUTPUT1<-list(NCASES1, NVALID1, NVALID2, NVALID12,
                  CORC, COR1$p.value,
                  mean(OUTPUT0), sd(OUTPUT0),
                  min(OUTPUT0), max(OUTPUT0),
                  mean(OUTPUT0<=CORC), mean(OUTPUT0>=CORC),
                  min(1,mean(OUTPUT0<=(-abs(CORC))) + mean(OUTPUT0>=(abs(CORC)))),
                  NPERM)
    names(OUTPUT1)<-c("NumberOfCases", "ValidCases_VEC1", "ValidCases_VEC2",
                      "ValidCases_Both",
                      "CorrCoeff","ClassicSign_2tailed",
                      "Mean_Perm","SD_Perm",
                      "Min_Perm","Max_Perm",
                      "AsSmall","AsLarge","AsExtreme",
                      "NumberOfPermutations")

    hist(OUTPUT0)
  }
  else
  {
    OUTPUT1<-list(NCASES1, NVALID1, NVALID2, NVALID12,
                  CORC, COR1$p.value)
    names(OUTPUT1)<-c("NumberOfCases", "NumberOfValidVEC1", "NumberOfValidVEC2", "NumberOfValidBoth",
                      "CorrelationCoef","SignClassic")
  }
  OUTPUT2<-t(t(unlist(OUTPUT1)))
  colnames(OUTPUT2)<-"Value"
  print(OUTPUT2)
  invisible(OUTPUT1)
}

#' Performs a regression analysis using a classic (and a permutation-based) significance test on a vector variable
#'
#' This function takes a vectors as a dependent variable and regresses it on a set of independent variables.
#' It returns the regression coefficients and the classic (two-tailed) significance tests.
#' If NPERM (number of permutations) is set to a value higher than 0 (default is 1000), permutations-based significance is also provided,
#' which will provide information about the regression coefficients for the permuted variables. A higher value than 1000 is generally advisable. Both the one-tailed ("AsSmall" or "AsLarge")
#' and the two-tailed ("AsExtreme") test are given. Note that cases with missing values are included in the permutation and this can impact the degrees of freedom for the regression analysis obtained from the permutations.
#' For example, suppose the same 20% of cases are missing for both variables (i.e., unit non-response), then permuting these missings will mean that a lot more cases will have a missing for at least one of the variables.
#' - Last updated: 23 April 2022.
#'
#' @param VEC1 A vector of values for the dependent variable.
#' @param DAT2 A data.frame containing the independent variables.
#' @param NPERM A positive integer value indicating the number of permutations to be used for the statistical test. If NPERM is set to NULL, NA or 0, no permutation-based significance test is performed.
#'
#' @return Returns an object containing information about the model and the regression coefficients.
#' @importFrom stats lm pf .lm.fit cor
#' @export

#' @author Filip Agneessens <filipagneessens2@gmail.com>
#'
#' @examples
#' VEC1<-c(2,4,2,3,3,1)
#' X1<-c(4,2,3,2,1,5)
#' X2<-c(1,0,0,1,1,0)
#' DAT2<-data.frame(X1,X2)
#' RESULTS1<-xRegression(VEC1,DAT2,NPERM=1000)
#' RESULTS1
#' xRegression(VEC1,DAT2,NPERM=0)

xRegression<-function(VEC1,DAT2,NPERM=1000)
{
  options(scipen=999)
  if(is.null(NPERM)) {NPERM<-0} # changed to 0 because NULL is not a value
  if(is.na(NPERM)) {NPERM<-0} # changed to 0 because NA is not a value

  if(!is.vector(VEC1))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "VEC1" is not a vector.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if(!is.data.frame(DAT2))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "DAT2" is not a data.frame.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if((ceiling(NPERM)!=floor(NPERM)) | (NPERM<0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The number of permutations ("NPERM") requested is not a positive integer value.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  NCASES1<-length(VEC1)
  DIM2<-dim(DAT2)

  if(NCASES1!=DIM2[1])
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The length for VEC1 is not the same as the number of rows for DAT2.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  REGR1<-lm(VEC1~.,data=DAT2)

  #MODEL INFO
  ADJRSQ<-summary(REGR1)$adj.r.squared
  RSQ<-summary(REGR1)$r.squared
  FTEST<-summary(REGR1)$fstatistic
  names(FTEST)<-NULL

  PMODEL<-pf(FTEST[1],FTEST[2],FTEST[3],lower.tail=F)

  NVALID1<-sum(!is.na(VEC1))
  NVALID2<-colSums(!is.na(DAT2))
  NVALID12<-sum(!is.na(rowSums(DAT2)*VEC1))
  CNAMES<-colnames(DAT2)

  #MODEL INFO
  OUTPUT1<-list()
  OUTPUT1$Model<-list(NCASES1, NVALID1, NVALID2, NVALID12, RSQ, ADJRSQ, FTEST[1], FTEST[2], FTEST[3], PMODEL, NPERM)
  names(OUTPUT1$Model)<-c("NumberOfCases", "NValidCases_DV", "NValidCases_IV",
                          "NValidCasesAll", "RSquared", "AdjustedRSquared", "FValue", "DF1", "DF2", "Significance", "NumberOfPermutations")

  #COEFF INFO
  Regr_Coeff<-summary(REGR1)$coefficients
  REGRC<-REGR1$coefficients

  # Only perform a permutation test if NPERM>0 (and not NULL or NA).
  if(NPERM>0)
  {
    #IV
    MAT2<-cbind(rep(1,length(VEC1)),as.matrix(DAT2))
    colnames(MAT2)<-c("Intercept",CNAMES)

    #PERMUTED REGRESSIONS
    OUTPUT0<-replicate(NPERM,.lm.fit(MAT2,sample(VEC1))$coefficients)

    #GET OUTPUT
    # FIRST STANDARD REGRESSION
    Mean_Perm<-rowMeans(OUTPUT0)
    SD_Perm<-apply(OUTPUT0,1,sd,na.rm = TRUE)
    Prop_AsSmall<-rowMeans(OUTPUT0<=REGRC)
    Prop_AsLarge<-rowMeans(OUTPUT0>=REGRC)
    Prop_AsExtreme<-rowMeans(OUTPUT0<=-abs(REGRC))+rowMeans(OUTPUT0<=abs(REGRC))

    OUTPUT1$Variables<-data.frame(Regr_Coeff, Mean_Perm, SD_Perm, Prop_AsSmall, Prop_AsLarge, Prop_AsExtreme)
    colnames(OUTPUT1$Variables)<-c("Coeff","StErr","t_value","Classic_Sign","Mean_Perm",
                                   "SD_Perm", "AsSmall", "AsLarge", "AsExtreme")
  }
  else
  {
    OUTPUT1$Variables<-data.frame(Regr_Coeff)
    colnames(OUTPUT1$Variables)<-c("Coeff","StErr","t_value","Classic_Sign")
  }

  cat('$Model\n')
  MOUT<-round(t(t(unlist(OUTPUT1$Model))),digits=5)
  colnames(MOUT)<-"Value"
  print(MOUT)
  cat('\n$Variables\n')
  print(round(OUTPUT1$Variables,digits=5))
  invisible(OUTPUT1)
  cat('\nClassic significance value based on assumption of independence of cases.\n')
}

#' Takes a one-mode network and randomly permutes both rows and columns
#'
#' This function takes a single one-mode network and randomly reorders (permutes) both rows and columns. This function is used for QAP-based analysis.
#' - Last updated: 23 April 2022.
#'
#' @param NET0 A network (currently only adjacency matrices are supported).
#'
#' @return A permuted adjacency matrix, where rows and columns are reodered in the same way.
#' @export
#' @author Filip Agneessens <filipagneessens2@gmail.com>
#'
#' @examples
#' MAT1<-matrix(c(2,4,2,3,3,1,3,4,1),3,3)
#' xPermuteQAP(MAT1)

xPermuteQAP<-function(NET0)
{
  PERMV<-sample(1:nrow(NET0))
  NET0[PERMV,PERMV]
}

#' Performs a correlation test using a classic (and a permutation-based) significance test between two network variables
#'
#' This function takes two networks among the same nodes and returns the correlation coefficient("CorrCoef") and the classic (two-tailed) significance test ("ClassicP_2T").
#' If NPERM (number of permutations) is set to a value higher than 0 (default is 1000), a quadratic assignment permutations-based (QAP) significance is also provided,
#' which will provide information about the correlations coefficients for the permuted networks (including the histogram). A higher value than 1000 is generally advisable. Both the one-tailed ("Prop_SmallerThan" or "Prop_GreaterThan")
#' and the two-tailed ("Prop_MoreExtreme") test are given. Note that cases with missing values are included in the permutation and this can impact the degrees of freedom for the correlations obtained for the permutations.
#' For example, suppose the same 20% of cases are missing for both variables (i.e., unit non-response), then permuting these missings will mean that a lot more cases will have a missing for at least one of both variables.
#' - Last updated: 23 April 2022.
#'
#' @param NET1 A first network.
#' @param NET2 A second network, among the same nodes as the first network.
#' @param NPERM A positive integer value indicating the number of permutations to be used for the statistical test. If NPERM is set to NULL, NA or 0, no permutation-based significance test is performed.
#' @param Directed Whether the analysis should consider ties in both directions (node i to j and j to i) as separate cases. Default is TRUE.
#' @param Loops Whether to allow self-loops.
#'
#' @return Prints a matrix with at least the following values: "NumberOfCases", "NumberOfValid1", "NumberOfValid2", "NumberOfValidBoth", "CorrCoef", "ClassicP_2T", and in case a permutation test is performed, in addition: Mean, StandardDev, Min_Perm, Max_Perm, Prop_SmallerThan, Prop_GreaterThan, Prop_MoreExtreme, Prop_SmallerThanABS, Prop_GreaterThanABS, NumberOfPermut. When using "<-" creates a list object with these respective values as separate elements.
#' @importFrom stats cor.test cor
#' @importFrom graphics hist
#' @export
#' @author Filip Agneessens <filipagneessens2@gmail.com>
#'
#' @examples
#' X1<-matrix(c(2,4,2,3,3,1,3,4,1),3,3)
#' X2<-matrix(c(4,1,5,2,1,4,2,2,3),3,3)
#' RESULS1<-xQAPCorrelation(X1,X2, NPERM=5000)
#' RESULS1
#' xQAPCorrelation(X1,X2,NPERM=0)

xQAPCorrelation<-function(NET1,NET2,NPERM=1000, Directed=TRUE, Loops=FALSE)
{
  options(scipen=999)
  if(is.null(NPERM)) {NPERM<-0} # changed to 0 because NULL is not a value
  if(is.na(NPERM)) {NPERM<-0} # changed to 0 because NA is not a value

  if(!is.matrix(NET1))

  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "NET1" is not a matrix.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if(!is.matrix(NET2))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "NET2" is not a matrix.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if((ceiling(NPERM)!=floor(NPERM)) | (NPERM<0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The number of permutations ("NPERM") requested is not a positive integer value.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  DIM1<-dim(NET1)
  DIM2<-dim(NET2)

  if(DIM1[1]!=DIM1[2])
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET1 is not a one-mode network.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(DIM2[1]!=DIM2[2])
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET2 is not a one-mode network.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(DIM1[1]!=DIM2[1])
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET1 and NET2 have a different number of nodes.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  vNET1<-NET1
  is.na(vNET1)<-0
  vNET2<-NET2
  is.na(vNET2)<-0
  if(Directed==FALSE & (sum(vNET1!=t(vNET1))>0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET1 is supposed to be undirected (as Directed=FALSE), yet the data in NET1 is not symmetrical.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(Directed==FALSE & (sum(vNET2!=t(vNET2))>0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET2 is supposed to be undirected (as Directed=FALSE), yet the data in NET2 is not symmetrical.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  NN1<-DIM1[1]
  NCASES<-NN1*(NN1-1)/(2-Directed)+Loops*NN1
  if(Loops==FALSE)
  {
    diag(NET1)<-NA
    diag(NET2)<-NA
  }
  NET1x<-NET1
  NET2x<-NET2

  if(Directed==FALSE)
  {
    NET1x[upper.tri(NET1x)]<-NA
    NET2x[upper.tri(NET2x)]<-NA
  }

  NVALID1<-sum(!is.na(NET1x))
  NVALID2<-sum(!is.na(NET2x))
  NVALID12<-sum(!is.na(NET1x*NET2x))

  COR1<-cor.test(c(NET1x),c(NET2x))
  CORC<-COR1$estimate
  names(CORC)<-NULL

  # only perform QAP analysis if the value for NPERM is not NULL and not NA and not 0.
  if(NPERM!=0)
  {
    OUTPUT0<-replicate(NPERM, cor(c(xPermuteQAP(NET1)),c(NET2), use = "pairwise.complete.obs"))
    OUTPUT1<-list(NN1, NCASES, NVALID1, NVALID2, NVALID12,
                  CORC, COR1$p.value,
                  mean(OUTPUT0), sd(OUTPUT0),
                  min(OUTPUT0), max(OUTPUT0),
                  mean(OUTPUT0<=CORC), mean(OUTPUT0>=CORC),
                  min(1,mean(OUTPUT0<=(-abs(CORC))) + mean(OUTPUT0>=(abs(CORC)))),
                  NPERM)
    names(OUTPUT1)<-c("NumberOfNodes", "NumberOfCases", "ValidCases_NET1", "ValidCases_NET2",
                      "ValidCases_Both",
                      "CorrCoef","ClassicSign_2tailed",
                      "Mean_Perm","SD_Perm",
                      "Min_Perm","Max_Perm",
                      "AsSmall","AsLarge","AsExtreme",
                      "NumberOfPermutations")
    hist(OUTPUT0)
  }
  else
  {
    OUTPUT1<-list(NN1, NCASES, NVALID1, NVALID2, NVALID12,
                  CORC, COR1$p.value)

    names(OUTPUT1)<-c("NumberOfNodes", "NumberOfCases", "ValidCases_NET1", "ValidCases_NET2",
                      "ValidCases_Both",
                      "CorrCoef","ClassicSign_2tailed")
  }

  OUTPUT2<-t(t(unlist(OUTPUT1)))
  colnames(OUTPUT2)<-"Value"
  print(OUTPUT2)
  invisible(OUTPUT1)
}


#' Performs a regression analysis using a classic (and a permutation-based) QAP significance test on a network variable
#'
#' This function takes a network as a dependent variable and regresses it on a set of independent networks. The approach uses each tie as a case.
#' It returns the regression coefficients and the classic (two-tailed) significance tests. Because of the independence of cases (ties) the classic significance test should not be relied on.
#' If NPERM (number of permutations) is set to a value higher than 0 (default is 1000), a quadratic assignment permutations-based (QAP) significance is also produced,
#' which will provide information about the regression coefficients for the permuted variables. A higher value than 1000 is generally advisable.
#' Both the one-tailed ("AsSmall" or "AsLarge") and the two-tailed ("AsExtreme") test are given. This function relies on the netlm function from the sna package.
#' - Last updated: 26 April 2022.
#'
#' @param NET1 A first network (adjacency matrix), capturing the values for the dependent variable.
#' @param LIST2 A list object containing the networks (adjacency matrices) which will be the independent variables.
#' @param NPERM A positive integer value indicating the number of permutations to be used for the statistical test. If NPERM is set to NULL, NA or 0, no permutation-based significance test is performed.
#' @param Directed Whether the analysis should consider ties in both directions (node i to j and j to i) as separate cases. Default is TRUE.
#' @param Method The type of permutation being used. The options are the Double Dekker semi-partialling ("qapspp"), permuting the dependent variable (Y) ("qapy") and independently permuting each of the independent variables (X) ("qapallx"). See the argument "nullhyp" in the netlm function for the sna packages for details.
#' @param Loops Whether to allow self-loops.
#'
#' @return Prints the regression output.
#' @importFrom sna netlm
#' @importFrom stats lm pf cor
#' @export
#'
#' @references
#' Chapter 14. Borgatti, S. P., Everett, M. G., Johnson, J. C., & Agneessens, F. (2022). Analyzing Social Networks Using R. SAGE.
#' Dekker, D., Krackhardt, D., and Snijders, T.A.B. (2007) Sensitivity of MRQAP tests to collinearity and autocorrelation conditions. Psychometrika 72, 563-581.
#' Krackhardt, D. (1988) Predicting with networks: Nonparametric multiple-regression analysis of dyadic data. Social Networks 10, 359-381.
#' Butts, Carter T. (2016). sna: Tools for Social Network Analysis. R package version 2.4.
#' Butts, Carter T. (2008). Social Network Analysis with sna. Journal of Statistical Software, 24(6).
#'
#' @seealso [xUCINET::xQAPCorrelation()]
#'
#' @examples
#' MATRIX1<-matrix(c(0,4,2,3, 3,0,4,1, 4,8,0,2, 6,2,6,0),4,4)
#' MATRIX2<-matrix(c(0,1,2,1, 3,0,3,7, 4,7,0,2, 5,2,2,0),4,4)
#' MATRIX3<-matrix(c(0,4,5,3, 2,0,2,8, 4,4,0,2, 3,5,3,0),4,4)
#'
#' LIST2<-list(MATRIX2,MATRIX3)
#' RESULTS1<-xQAPRegression(MATRIX1,LIST2,NPERM=1000)
#' RESULTS1

xQAPRegression<-function(NET1, LIST2, NPERM=1000, Directed=TRUE, Method="qapspp", Loops=FALSE)
{
  options(scipen=999)

  if(is.null(NPERM)) {NPERM<-0} # changed to 0 because NULL is not a value
  if(is.na(NPERM)) {NPERM<-0} # changed to 0 because NA is not a value
  if(!Method %in% c("qapspp", "qapy", "qapallx"))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "Method" should be one of the following: "qapspp", "qapy", "qapallx".\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(!is.matrix(NET1))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "NET1" is not a matrix.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if(!is.list(LIST2))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "LIST2" is not a list.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  if((ceiling(NPERM)!=floor(NPERM)) | (NPERM<0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The number of permutations ("NPERM") requested is not a positive integer value.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  DIM1<-dim(NET1)
  LENGTH2<-length(LIST2)

  #check if any missings in NET1
  if(sum(is.na(NET1))!=0)
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET1 contains missings, but this function does not currently support missing values.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  #check if any missings in LIST2
  for(Nk in 1:LENGTH2)
  {
    if(sum(is.na(LIST2[[Nk]]))!=0)
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: Element ',Nk,' of LIST2 contains missings, but this function does not currently support missing values.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  #check if NET1 is square (one-mode) network
  if(DIM1[1]!=DIM1[2])
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET1 is not a one-mode network.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  #Matrix class and dimensions
  for(Nk in 1:LENGTH2)
  {

    if(!is.matrix(LIST2[[Nk]]))
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: Element ',Nk,' of "LIST2" is not a matrix.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
    DIMk<-dim(LIST2[[Nk]])

    if(DIMk[1]!=DIMk[2])
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: Element ',Nk,' of "LIST2" is not a one-mode network.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
    if(DIM1[1]!=DIMk[1])
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The number of rows in the ',Nk,'th element of "LIST2" is not equal to the number of rows in the DV "NET1".\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
    if(DIM1[2]!=DIMk[2])
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The number of columns in the ',Nk,'th element of "LIST2" is not equal to the number of columns in the DV "NET1".\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  #Directed?
  if(Directed==FALSE & (sum(NET1!=t(NET1))>0))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: NET1 is supposed to be undirected (as Directed=FALSE), yet the data in NET1 is not symmetrical.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  for(Nk in 1:LENGTH2)
  {
    if(Directed==FALSE & (sum(LIST2[[Nk]]!=t(LIST2[[Nk]]))>0))
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: Element ', Nk ,' of LIST2 is supposed to be undirected (as Directed=FALSE), yet the data in this matrix is not symmetrical.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  NN1<-DIM1[1]
  NCASES<-NN1*(NN1-1)/(2-Directed)+Loops*NN1
  NET1x<-NET1
  LIST2x<-LIST2

  if(Loops==FALSE)
  {
    diag(NET1x)<-NA
    for (Nk in 1:length(LIST2))
    {
      diag(LIST2x[[Nk]])<-NA
    }
  }

  if(Directed==FALSE)
  {
    NET1x[upper.tri(NET1x)]<-NA
    for (Nk in 1:length(LIST2))
    {
      LIST2x[[Nk]][upper.tri(LIST2x[[Nk]])]<-NA
    }
  }

  NET1x<-as.vector(NET1x)
  MAT2x<-data.frame(matrix(unlist(LIST2x),NN1*NN1,length(LIST2)))
  REGR1<-lm(NET1x~.,data=MAT2x)

  NVALID1<-sum(!is.na(NET1x))  # exclude diagonal if needed
  NVALID2<-colSums(!is.na(MAT2x))   # exclude diagonal if needed
  NVALID12<-sum(!is.na(rowSums(MAT2x)*NET1x))   # exclude diagonal if needed

  Regr_Coeff<-summary(REGR1)$coefficients
  Adjusted_R_Squared<-summary(REGR1)$adj.r.squared
  R_squared<-summary(REGR1)$r.squared
  F_test<-summary(REGR1)$fstatistic
  names(F_test)<-NULL
  Sign_Model<-pf(F_test[1],F_test[2],F_test[3],lower.tail=F)


  OUTPUT1<-list()
  OUTPUT1$Model<-list(NN1, NCASES, NVALID1, NVALID2, NVALID12, R_squared, Adjusted_R_Squared, F_test[1], F_test[2], F_test[3], Sign_Model, NPERM)
  names(OUTPUT1$Model)<-c("NumberOfNodes", "NumberOfCases", "NValidCases_DV", "NValidCases_IV",
                          "NValidCasesAll", "RSquared", "AdjustedRSquared", "FValue", "DF1", "DF2", "Significance", "NumberOfPermutations")

  #Standardized coefficients
  COEF1<-summary(REGR1)$coef[-1,1]
  SDX<-apply(REGR1$model[-1],2,sd,na.rm = TRUE)
  SDY<-sd(REGR1$model[[1]])
  BETA<-c(0,COEF1*SDX/SDY)

  # only perform QAP analysis if the value for NPERM is not NULL and not NA and not 0.
  if(NPERM!=0)
  {
    #Now use netlm in sna, using the correct specifications
    if(Directed==TRUE)
    {
      mode.sna<-"digraph"
    }
    else
    {
      mode.sna<-"graph"
    }

    REGR_P<-sna::netlm(NET1, LIST2, mode=mode.sna, nullhyp=Method, reps=NPERM, diag=Loops)
    OUTPUT1$MP<-summary(REGR_P)

    TST<-REGR_P$tstat
    OUTPUT0<-t(REGR_P$dist)
    Mean_Perm<-rowMeans(OUTPUT0)
    SD_Perm<-apply(OUTPUT0,1,sd,na.rm = TRUE)
    Min_Perm<-apply(OUTPUT0,1,min,na.rm = TRUE)
    Max_Perm<-apply(OUTPUT0,1,max,na.rm = TRUE)
    Prop_AsSmall<-rowMeans(OUTPUT0<=TST)
    Prop_AsLarge<-rowMeans(OUTPUT0>=TST)
    Prop_AsSmallABS<-rowMeans(OUTPUT0<=(-abs(TST)))
    Prop_AsLargeABS<-rowMeans(OUTPUT0>=(abs(TST)))
    Prop_AsExtreme<-Prop_AsSmallABS+Prop_AsLargeABS

    OUTPUT1$Variables<-data.frame(BETA, Regr_Coeff, Mean_Perm, SD_Perm, Prop_AsSmall, Prop_AsLarge, Prop_AsExtreme)
    colnames(OUTPUT1$Variables)<-c("Stand_Coeff","Coeff","StErr","t_value","[Sign*]","Mean_Perm",
                                   "SD_Perm", "AsSmall", "AsLarge", "AsExtreme")
  }
  else
  {
    OUTPUT1$Variables<-data.frame(BETA,Regr_Coeff)
    colnames(OUTPUT1$Variables)<-c("Stand_Coeff","Coeff","StErr","t_value","[Sign*]")
  }
  cat('\n$Model\n')
  MOUT<-round(t(t(unlist(OUTPUT1$Model))),digits=5)
  colnames(MOUT)<-"Value"
  print(MOUT)
  cat('\n$Variables\n')
  print(round(OUTPUT1$Variables,digits=5))
  cat('\nSign* = classic significance value based on assumption of independence of cases. Not to be used.\n')
  invisible(OUTPUT1)
}
