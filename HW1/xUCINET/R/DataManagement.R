#' Creating a "dataproject" containing network data
#'
#' This function creates a new "dataproject" list object, and adds the network (NETFILE1) and extra information about the network.
#' - Last updated: 28 April 2022.
#'
#' @param GeneralDescription Detailed information about the project as a whole.
#' @param NetworkName The name to be given to the network.
#' @param NETFILE1 A single social network file. This could be a csv file or an R object. The input type is defined in "FileType", while the format of the input data is defined in "InFormatType".
#' @param FileType The type of data that is being imported. Default is "csv". Other options are "Robject" when the data is already an object in R, and "csv2". Note that by default a csv file use a point (".") for the decimal point and a comma for the separator in most countries. However, in some Western European countries a comma is used for the decimal point and a semicolon for the separator, which will require "csv2".
#' @param InFormatType The format used to store the network data in a csv file. Currently only reads in an adjacency matrix "AdjMat".
#' @param NetworkDescription Detailed information about the network being imported.
#' @param Mode Character (string), indicating a name for the mode of the network. For two-mode networks this would need to contain the names for both modes. If it is a one-mode network, and Mode=NA then the mode name "A" will be given.
#' @param Directed Logical, whether the network should be considered directed or not. The default is TRUE.
#' @param Loops Logical, whether the network's self-nominations (loops) should be considered. By default this is FALSE and loops are not considered.
#' @param Values Character string, indicating what type of values the network contains (e.g., "Binary","Nominal","Ordinal","Rank","Interval","Ratio").
#' @param Class The way the data is being stored. Currently the network will be stored as an adjacency matrix in an object of class "matrix".
#' @param References Any references.
#'
#' @return A "dataproject" object containing network data.
#' @importFrom utils read.csv read.csv2
#' @export
#'
#' @author Filip Agneessens, \email{filipagneessens2@@gmail.com}
#'
#' @examples
#' ## Create a network with ordinal categories (0='no relation', 1='acquaintance', 2='friend',
#' ##    3='good friend') as an (adjacency) matrix:
#' Friendship3<-matrix(c(0,1,1, 2,0,3, 0,2,0),3,3, byrow=TRUE)
#' ## Add node names to the matrix:
#' rownames(Friendship3)<-c("a","b","c")
#' colnames(Friendship3)<-rownames(Friendship3)
#'
#' ## Now create the "dataproject", which we call "SchoolClass_Project1"
#' SchoolClass_Project1<-xCreateProject(GeneralDescription="Dataset among 3 students",
#'      NetworkName="Friendship",
#'      NETFILE1=Friendship3,
#'      FileType="Robject",
#'      InFormatType="AdjMat",
#'      NetworkDescription="Strength (0=no relation, 1=acquaintance, 2=friend, 3=good friend)",
#'      Mode=c("Students"),
#'      Directed=TRUE,
#'      Loops=FALSE,
#'      Values="Ordinal",
#'      Class="matrix",
#'      References="No references")
#' ## Check the newly created "dataproject" object
#' SchoolClass_Project1
#'
xCreateProject<-function(GeneralDescription=NULL,
                         NetworkName=NULL,
                         NETFILE1,
                         FileType=NULL,
                         InFormatType="AdjMat",
                         NetworkDescription=NULL,
                         Mode=NA,
                         Directed=TRUE,
                         Loops=FALSE,
                         Values=NA,
                         Class="matrix",
                         References)
{
  cat('\n ------ FUNCTION: xCreateProject ---------------------------\n\n')
  if(is.null(NetworkName))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: A name for the network (NetworkName) should be specified.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(NetworkDescription))
  {
    NetworkDescription<-"No information about this network."
  }

  if(is.null(FileType))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The FileType should be "csv", "csv2", or "Robject"\n',
         '##        Choose among: [FileType="csv"], [FileType="csv2"], and [FileType="Robject"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if((FileType!="csv") & (FileType!="csv2") & (FileType!="Robject"))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The FileType should be "csv", "csv2", or "Robject"\n',
         '##        Choose among: [FileType="csv"], [FileType="csv2"], and [FileType="Robject"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(InFormatType))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The InFormatType should be "AdjMat"',
         '##        Choose: [InFormatType="AdjMat"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(InFormatType!="AdjMat")
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Currently the InFormatType should be "AdjMat"',
         '##        Choose: [InFormatType="AdjMat"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(Class))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Currently this function only supports importing networks into a class "matrix"\n',
         '##        Choose: [Class="matrix"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(Class!=c("matrix"))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Currently this function only supports importing networks into a class "matrix"\n',
         '##        Choose: [Class="matrix"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(Mode))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Mode should be a vector of length 1 (e.g.,"Person"), or 2 (e.g.,c("Person","Event"))\n',
         '##        Adjust the length of Mode, e.g. [Mode="A"], or [Mode=c("A","B")].\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(length(Mode)>2)
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Mode should be a vector of length 1 (e.g.,"Person"), or 2 (e.g.,c("Person","Event"))\n',
         '##        Adjust the length of Mode, e.g. [Mode="A"], or [Mode=c("A","B")].\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(length(Mode)==2)
  {
    if(Mode[1]==Mode[2])
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The two elements of the vector provided for Mode are the same.\n',
           '##        Provide unique names for both, such as [Mode=c("Person","Event")] if it is a two-mode network.\n',
           '##        or provide a single name, such as [Mode=c("Person")] if it is a one-mode network.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  cat('\n - Basic checks performed on the argument for "xCreateProject". \n')

  # Importing a csv file or object
  if(FileType=="csv")
  {
    #check for errors
    NET1<-as.matrix(read.csv(NETFILE1, stringsAsFactors=FALSE, row.names=1))
  }

  if(FileType=="csv2")
  {
    #check for errors
    NET1<-as.matrix(read.csv2(NETFILE1, stringsAsFactors=FALSE, row.names=1))
  }

  if(FileType=="Robject")
  {
    #check for errors
    NET1<-NETFILE1
  }

  cat('\n - Data imported: [', NETFILE1, '] and named as: [', NetworkName,']\n')

  # Add row and column names if absent
  if((is.null(rownames(NET1))) & (is.null(colnames(NET1))))
  {
    if(nrow(NET1)!=ncol(NET1))
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      CNames<-c(1:ncol(NET1))
      colnames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for either mode found. We will use "a1", "a2", etc. for the first mode and "e1", "e2", etc. for the second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==2)
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      CNames<-c(1:ncol(NET1))
      colnames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for either mode found. We will use "a1", "a2", etc. for the first mode and "e1", "e2", etc. for the second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==1)
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      colnames(NET1)<-rownames(NET1)
      cat('\n No names for the nodes found. We will use "a1", "a2", etc.\n')
    }
  }

  # Add row names absent
  if((is.null(rownames(NET1))) & (!is.null(colnames(NET1))))
  {
    if(nrow(NET1)!=ncol(NET1))
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      cat('\n No names for rows mode found. We will use "a1", "a2", etc. for the first mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==2)
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      cat('\n No names for rows mode found. We will use "a1", "a2", etc. for the first mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==1)
    {
      rownames(NET1)<-colnames(NET1)
      cat('\n No names for the rows found. We copied the column names.\n')
    }
  }

  # Add column names absent
  if((!is.null(rownames(NET1))) & (is.null(colnames(NET1))))
  {
    if(nrow(NET1)!=ncol(NET1))
    {
      CNames<-c(1:ncol(NET1))
      rownames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for column mode found. We will use "e1", "e2", etc. for this second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==2)
    {
      CNames<-c(1:ncol(NET1))
      rownames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for column mode found. We will use "e1", "e2", etc. for this second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==1)
    {
      colnames(NET1)<-rownames(NET1)
      cat('\n No names for the columns found. We copied the row names.\n')
    }
  }

  #-----
  if(length(Mode)==1)
  {
    #if Mode==NA
    if(is.na(Mode))
    {
      if(nrow(NET1)!=ncol(NET1))
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: It seems we are dealing with a two-mode network as the number of rows and columns\n',
             '##        are not the same. However, the default option of Mode=NA is provided as input.\n',
             '##        If this is a two-mode network, provide unique names for both modes, such as [Mode=c("A","B")].\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }

      if(nrow(NET1)==ncol(NET1))
      {
        #if Mode==NA and nrow==ncol but row/colnames are different
        if(sum(rownames(NET1)!=colnames(NET1))>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: It is unclear whether we are dealing with a one-mode or a two-mode network.\n',
               '##        The row and column names for the network do not correspond,\n',
               '##        however, the number of rows and columns are the same and "Mode=NA".\n',
               '##        Check row/column names and/or specify "Mode".\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
        else {Mode<-c("A")}
      }
    }
  }

  if (length(Mode)==1)
  {
    #this excludes that Mode is still NA, which has been dealt with before and replaced.
    if(nrow(NET1)!=ncol(NET1))
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The argument "Mode" contains a vector of length 1.\n',
           '##        However the network is not one-mode (not square).\n',
           '##        If it is a two-mode network, specify both modes using the format Mode=c("A","B")',
           '##        If it is a one-mode network, check the input data for mistakes.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
    if(sum(rownames(NET1)!=colnames(NET1))>0)
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The row and column names are not equal, yet the network is supposed to be one-mode.\n',
           '##        Check row/column names and/or adjust "Mode" if this is a two-mode network.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  #-Directed yet Directed=FALSE---
  if(length(Mode)==1)
  {
    if(Directed==FALSE)
    {
      if(sum(NET1!=t(NET1))>0){
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: Directed is specified as FALSE, yet the network is not symmetrical.\n',
             '##        Check the network or set "Directed=TRUE" if it is a directed network.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }
    }
  }

  # Create a dataframe that only contains the names of nodes as an attribute file

  if(length(Mode)==1)
  {
    Mode1<-Mode
    Mode2<-Mode
    ATTRIBUTEFILE1<-data.frame(rownames(NET1))
    colnames(ATTRIBUTEFILE1)<-"NodeName"
  }
  if(length(Mode)==2)
  {
    Mode1<-Mode[1]
    Mode2<-Mode[2]
    ATTRIBUTEFILE1<-list(data.frame(rownames(NET1)),data.frame(colnames(NET1)))
    names(ATTRIBUTEFILE1)<-c(Mode1,Mode2)
    colnames(ATTRIBUTEFILE1[[1]])<-"NodeName"
    colnames(ATTRIBUTEFILE1[[2]])<-"NodeName"
  }

  # Creating "dataproject" object
  AttributesDescription<-data.frame("NodeName",Mode1,paste("Names of the nodes for mode",Mode1))
  if(length(Mode)==2)
  {
    AttributesDescription<-rbind(AttributesDescription,c("NodeName",Mode2,paste("Names of the nodes for mode",Mode2)))
  }
  colnames(AttributesDescription)<-c("Variable","Mode","Details")
  NetworksDescription<-data.frame(NetworkName,NetworkDescription)
  colnames(NetworksDescription)<-c("NetworkName","Details")
  OUTPUT1<-list(GeneralDescription,Mode,AttributesDescription,NetworksDescription,References)
  names(OUTPUT1)<-c("GeneralDescription","Modes","AttributesDescription","NetworksDescription","References")
  OUTPUT2<-data.frame(NetworkName,Mode1,Mode2,Directed,Loops,Values,Class)
  names(OUTPUT2)<-c("NetworkName","ModeSender","ModeReceiver","Directed","Loops","Values","Class")

  PROJECT1<-list(OUTPUT1,ATTRIBUTEFILE1,OUTPUT2,NET1)
  names(PROJECT1)<-c("ProjectInfo","Attributes","NetworkInfo",NetworkName)
  # Output project
  cat(' ---------------------------------------------------------------\n\n')

  PROJECT1
}

#' Add an additional network to an existing "dataproject"
#'
#' This function adds a network and information about that network to an existing "dataproject".
#' - Last updated: 28 April 2022.
#'
#' @param ProjectName The name of the "dataproject", to which the network needs to be added.
#' @param NetworkName The name to be given to the new network.
#' @param NETFILE1 A single social network file. This could be a csv file or an R object. The input type is defined in "FileType", while the format of the input data is defined in "InFormatType".
#' @param FileType The type of data that is being imported. Default is "csv". Other options are "Robject" when the data is already an object in R, and "csv2". Note that by default a csv file use a point (".") for the decimal point and a comma for the separator in most countries. However, in some Western European countries a comma is used for the decimal point and a semicolon for the separator, which will require "csv2".
#' @param InFormatType The format used to store the network data in a csv file. Currently only reads in an adjacency matrix "AdjMat".
#' @param NetworkDescription Detailed information about the network being imported.
#' @param Mode Character (string), indicating a name for the mode of the network. For two-mode networks this would need to contain the names for both modes. If it is a one-mode network, and Mode=NA then the mode name "A" will be assumed.
#' @param Directed Logical, whether the network should be considered directed or not. The default is TRUE.
#' @param Loops Logical, whether the network's self-nominations (loops) should be considered. By default this is FALSE and loops are not considered.
#' @param Values Character string, indicating what type of values the network contains (e.g., "Binary","Nominal","Ordinal","Rank","Interval","Ratio").
#' @param Class The way the data is being stored. Currently the network will be stored as an adjacency matrix in an object of class "matrix".
#'
#' @return A "dataproject" object containing network data.
#' @importFrom utils read.csv read.csv2
#' @export
#'
#' @author Filip Agneessens, \email{filipagneessens2@@gmail.com}
#'
#' @examples
#' ## We can consider an existing network and add the transposed
#' ## Create a network with ordinal categories (0='no relation', 1='acquaintance', 2='friend',
#' ##    3='good friend') as an (adjacency) matrix:
#' Friendship3<-matrix(c(0,1,1, 2,0,3, 0,2,0),3,3, byrow=TRUE)
#' ## Add node names to the matrix:
#' rownames(Friendship3)<-c("a","b","c")
#' colnames(Friendship3)<-rownames(Friendship3)
#'
#' ## Now create the "dataproject", which we call "SchoolClass_Project1"
#' SchoolClass_Project1<-xCreateProject(GeneralDescription="Dataset among 3 students",
#'      NetworkName="Friendship",
#'      NETFILE1=Friendship3,
#'      FileType="Robject",
#'      InFormatType="AdjMat",
#'      NetworkDescription="Strength (0=no relation, 1=acquaintance, 2=friend, 3=good friend)",
#'      Mode=c("Students"),
#'      Directed=TRUE,
#'      Loops=FALSE,
#'      Values="Ordinal",
#'      Class="matrix",
#'      References="No references")
#' ## Check the newly created "dataproject" object
#' SchoolClass_Project1
#'
#' ## Now create the transposed
#' SchoolClass_Project1<-xAddToProject(ProjectName=SchoolClass_Project1,
#'         NetworkName="Friendship_Tr",
#'         NETFILE1=t(SchoolClass_Project1$Friendship),
#'         FileType="Robject",
#'         InFormatType="AdjMat",
#'         NetworkDescription="Transposed of the friendship network",
#'         Mode=c("Students"),
#'         Directed=TRUE,
#'         Loops=FALSE,
#'         Values="Ordinal",
#'         Class="matrix")
#' ## Check the new result:
#' SchoolClass_Project1

xAddToProject<-function(ProjectName,
                        NetworkName=NULL,
                        NETFILE1=NULL,
                        FileType=NULL,
                        InFormatType="AdjMat",
                        NetworkDescription=NULL,
                        Mode=NA,
                        Directed=TRUE,
                        Loops=FALSE,
                        Values=NA,
                        Class="matrix")
{
  cat('\n ------ FUNCTION: xAddToProject ---------------------------\n\n')
  ProjectName1<-deparse(substitute(ProjectName))

  # "dataproject"
  #check that "dataproject" with ProjectName exists
  if(class(try(get(ProjectName1), silent=TRUE))=="try-error")
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "ProjectName" not found. Make sure the name of the "dataproject" object exists,\n',
         '##        and does not contain typos, and that it is placed in between quotes "..." .\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  cat("ProjectName:",ProjectName1)

  #now (if it exists) get the "dataproject"
  PROJECT1<-get(ProjectName1)

  #length of list object is 3 more than the number of networks: [1] ProjectInfo, [2] Attributes, [3] NetworkInfo, [4] "Network1", ...
  NLength<-length(PROJECT1)
  cat("Number of existing networks:",(NLength - 3))

  if(is.null(NetworkName))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: A name for the network to be added (NetworkName) should be specified.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(NetworkName %in% names(PROJECT1)) {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The name for the new network', NetworkName, 'is the same as an existing one.\n',
         '##        Change the [NetworkName=', NetworkName,'] to a new name,\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(NetworkDescription))
  {
    NetworkDescription<-"No information about this network."
  }

  cat('\nNetworkName:',NetworkName)

  #Check if NETFILE1 is provided
  if(is.null(NETFILE1)) {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: No NETFILE1 provided\n',
         '##        Add the file name using [NETFILE1=...].\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(FileType))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The FileType should be "csv", "csv2", or "Robject"\n',
         '##        Choose among: [FileType="csv"], [FileType="csv2"], and [FileType="Robject"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if((FileType!="csv") & (FileType!="csv2") & (FileType!="Robject"))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The FileType should be "csv", "csv2", or "Robject"\n',
         '##        Choose among: [FileType="csv"], [FileType="csv2"], and [FileType="Robject"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(InFormatType))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The InFormatType should be "AdjMat"',
         '##        Choose: [InFormatType="AdjMat"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(InFormatType!="AdjMat")
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Currently the InFormatType should be "AdjMat"',
         '##        Choose: [InFormatType="AdjMat"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(Class))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Currently this function only supports importing networks into a class "matrix"\n',
         '##        Choose: [Class="matrix"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(Class!=c("matrix"))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Currently this function only supports importing networks into a class "matrix"\n',
         '##        Choose: [Class="matrix"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(is.null(Mode))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Mode should be a vector of length 1 (e.g.,"Person"), or 2 (e.g.,c("Person","Event"))\n',
         '##        Adjust the length of Mode, e.g. [Mode="A"], or [Mode=c("A","B")].\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(length(Mode)>2)
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Mode should be a vector of length 1 (e.g.,"Person"), or 2 (e.g.,c("Person","Event"))\n',
         '##        Adjust the length of Mode, e.g. [Mode="A"], or [Mode=c("A","B")].\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  if(length(Mode)==2)
  {
    if(Mode[1]==Mode[2])
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The two elements of the vector provided for Mode are the same.\n',
           '##        Provide unique names for both, such as [Mode=c("Person","Event")] if it is a two-mode network.\n',
           '##        or provide a single name, such as [Mode=c("Person")] if it is a one-mode network.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  cat('\n - Basic checks performed on the argument for "xAddToProject". \n')

  # Importing a csv file or object
  if(FileType=="csv")
  {
    #check for errors
    NET1<-as.matrix(read.csv(NETFILE1, stringsAsFactors=FALSE, row.names=1))
  }

  if(FileType=="csv2")
  {
    #check for errors
    NET1<-as.matrix(read.csv2(NETFILE1, stringsAsFactors=FALSE, row.names=1))
  }

  if(FileType=="Robject")
  {
    #check for errors
    NET1<-NETFILE1
  }

  cat('\n - Data imported: [', NETFILE1, '] and named: [', NetworkName,']\n')

  # Add row and column names if absent
  if((is.null(rownames(NET1))) & (is.null(colnames(NET1))))
  {
    if(nrow(NET1)!=ncol(NET1))
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      CNames<-c(1:ncol(NET1))
      colnames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for either mode found. We will use "a1", "a2", etc. for the first mode and "e1", "e2", etc. for the second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==2)
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      CNames<-c(1:ncol(NET1))
      colnames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for either mode found. We will use "a1", "a2", etc. for the first mode and "e1", "e2", etc. for the second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==1)
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      colnames(NET1)<-rownames(NET1)
      cat('\n No names for the nodes found. We will use "a1", "a2", etc.\n')
    }
  }

  # Add row names absent
  if((is.null(rownames(NET1))) & (!is.null(colnames(NET1))))
  {
    if(nrow(NET1)!=ncol(NET1))
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      cat('\n No names for rows mode found. We will use "a1", "a2", etc. for the first mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==2)
    {
      RNames<-c(1:nrow(NET1))
      rownames(NET1)<-paste("a",RNames, sep="")
      cat('\n No names for rows mode found. We will use "a1", "a2", etc. for the first mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==1)
    {
      rownames(NET1)<-colnames(NET1)
      cat('\n No names for the rows found. We copied the column names.\n')
    }
  }

  # Add column names absent
  if((!is.null(rownames(NET1))) & (is.null(colnames(NET1))))
  {
    if(nrow(NET1)!=ncol(NET1))
    {
      CNames<-c(1:ncol(NET1))
      rownames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for column mode found. We will use "e1", "e2", etc. for this second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==2)
    {
      CNames<-c(1:ncol(NET1))
      rownames(NET1)<-paste("e",CNames, sep="")
      cat('\n No names for column mode found. We will use "e1", "e2", etc. for this second mode.\n')
    }
    if((nrow(NET1)==ncol(NET1)) & length(Mode)==1)
    {
      colnames(NET1)<-rownames(NET1)
      cat('\n No names for the columns found. We copied the row names.\n')
    }
  }

  #-----
  if (length(Mode)==1)
  {
    #if Mode==NA
    if (is.na(Mode))
    {
      if(nrow(NET1)!=ncol(NET1))
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: It seems we are dealing with a two-mode network as the number of rows and columns\n',
             '##        are not the same. However, the default option of Mode=NA is provided as input.\n',
             '##        If this is a two-mode network, provide unique names for both modes, such as [Mode=c("A","B")].\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }

      if(nrow(NET1)==ncol(NET1))
      {
        #if Mode==NA and nrow==ncol but row/colnames are different
        if(sum(rownames(NET1)!=colnames(NET1))>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: It is unclear whether we are dealing with a one-mode or a two-mode network.\n',
               '##        The row and column names for the network do not correspond,\n',
               '##        however, the number of rows and columns are the same and Mode=NA.\n',
               '##        Check row/column names and/or specify "Mode".\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
        else {Mode<-c("A")}
      }
    }
  }

  if (length(Mode)==1)
  {
    #this excludes that is still Mode==NA, which has been dealt with before and replaced
    if(nrow(NET1)!=ncol(NET1))
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The argument "Mode" contains a vector of length 1.\n',
           '##        However the network is not one-mode (not square).\n',
           '##        If it is a two-mode network, specify both modes using the format Mode=c("A","B")',
           '##        If it is a one-mode network, check the input data for mistakes.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
    if(sum(rownames(NET1)!=colnames(NET1))>0)
    {
      stop('\n',
           '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
           '## ERROR: The row and column names are not equal, yet the network is supposed to be one-mode.\n',
           '##        Check row/column names and/or adjust "Mode" if this is a two-mode network.\n',
           '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
    }
  }

  #-Directed yet Directed=FALSE---
  if(length(Mode)==1)
  {
    if(Directed==FALSE)
    {
      if(sum(NET1!=t(NET1))>0){
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: Directed is specified as FALSE, yet the network is not symmetrical.\n',
             '##        Check the network or set "Directed=TRUE" if it is a directed network.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')

      }
    }
  }

  #---------------------------------------------------------------------------------------
  # Check names of new network against existing and add empty attribute file if needed.
  # If new network is one-mode
  if(length(Mode)==1)
  {
    Mode1<-Mode
    Mode2<-Mode
    if(length(PROJECT1$ProjectInfo$Modes)==1)
    {
      # If one existing mode -> if same as existing one-mode
      if(Mode==PROJECT1$ProjectInfo$Modes)
      {
        if(sum(rownames(NET1)!=PROJECT1$Attributes$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The row names of the network are not the same as those for the existing network.\n',
               '##        Check the names of the network or correctly specify the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
      }

      # If one existing mode -> if not the same as existing one-mode
      if(Mode!=PROJECT1$ProjectInfo$Modes)
      {
        # Make new mode-specific attribute
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode)
        # Create new dataframe
        ATTRIBUTEFILE1<-data.frame(rownames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        PROJECT1$Attributes<-list(PROJECT1$Attributes,ATTRIBUTEFILE1)
        names(PROJECT1$Attributes)<-PROJECT1$ProjectInfo$Modes
        #
        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode,paste("Names of the nodes for mode",Mode)))
      }
    }

    # If we already have more than one mode in the dataset
    if(length(PROJECT1$ProjectInfo$Modes)>1)
    {
      if(Mode %in% PROJECT1$ProjectInfo$Modes)
      {
        if(sum(rownames(NET1)!=PROJECT1$Attributes$Mode$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The row names of the network are not the same as those for the existing network.\n',
               '##        Check the names of the network or correctly specify the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
      }
      if(!(Mode %in% PROJECT1$ProjectInfo$Modes))
      {
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode)

        #create a new data.frame for the attributes for a new mode
        ATTRIBUTEFILE1<-data.frame(rownames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        PROJECT1$Attributes$NEW<-ATTRIBUTEFILE1
        names(PROJECT1$Attributes)[length(PROJECT1$Attributes)]<-Mode

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode,paste("Names of the nodes for mode",Mode)))
      }
    }
  }

  # New two-mode network
  if(length(Mode)==2)
  {
    Mode1<-Mode[1]
    Mode2<-Mode[2]

    if(length(PROJECT1$ProjectInfo$Modes)==1)
    {
      if(Mode1==PROJECT1$ProjectInfo$Modes)
      {
        if(sum(rownames(NET1)!=PROJECT1$Attributes$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The row names of the network are not the same as those for the existing network.\n',
               '##        Check the row names of the network or correctly specify the first value of the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
        #This means that Mode2 (column) is new
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode2)

        ATTRIBUTEFILE1<-data.frame(colnames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        PROJECT1$Attributes<-list(PROJECT1$Attributes,ATTRIBUTEFILE1)
        names(PROJECT1$Attributes)<-PROJECT1$ProjectInfo$Modes

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode2,paste("Names of the nodes for mode",Mode2)))
      }

      if(Mode2==PROJECT1$ProjectInfo$Modes)
      {
        if(sum(colnames(NET1)!=PROJECT1$Attributes$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The column names of the network are not the same as those for the existing network.\n',
               '##        Check the column names of the network or correctly specify the second value of the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
        #This means that Mode1 (rows) is new
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode1)

        ATTRIBUTEFILE1<-data.frame(rownames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        PROJECT1$Attributes<-list(PROJECT1$Attributes,ATTRIBUTEFILE1)
        names(PROJECT1$Attributes)<-PROJECT1$ProjectInfo$Modes

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode1,paste("Names of the nodes for mode",Mode1)))
      }

      if((Mode2!=PROJECT1$ProjectInfo$Modes) & (Mode1!=PROJECT1$ProjectInfo$Modes))
      {
        #This means that Mode1 and Mode2 are new
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode1,Mode2)

        ATTRIBUTEFILE1<-data.frame(rownames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        ATTRIBUTEFILE2<-data.frame(colnames(NET1))
        colnames(ATTRIBUTEFILE2)<-"NodeName"
        PROJECT1$Attributes<-list(PROJECT1$Attributes,ATTRIBUTEFILE1,ATTRIBUTEFILE2)
        names(PROJECT1$Attributes)<-PROJECT1$ProjectInfo$Modes

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode1,paste("Names of the nodes for mode",Mode1)))
        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode2,paste("Names of the nodes for mode",Mode2)))
      }
    }

    # If multiple modes already exist:
    if(length(PROJECT1$ProjectInfo$Modes)>1)
    {
      if((Mode1%in%PROJECT1$ProjectInfo$Modes) & (!Mode2%in%PROJECT1$ProjectInfo$Modes))
      {
        POSITION_Mode1<-sum((Mode1==names(PROJECT1$Attributes))*c(1:length(PROJECT1$Attributes)))

        if(sum(rownames(NET1)!=PROJECT1$Attributes[[POSITION_Mode1]]$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The row names of the network are not the same as those for the existing network.\n',
               '##        Check the row names of the network or correctly specify the first value of the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
        #This means that Mode2 (column) is new
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode2)

        ATTRIBUTEFILE1<-data.frame(colnames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        PROJECT1$Attributes$NEW<-ATTRIBUTEFILE1
        names(PROJECT1$Attributes)[length(PROJECT1$Attributes)]<-Mode2

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode2,paste("Names of the nodes for mode",Mode2)))
      }

      if((!Mode1%in%PROJECT1$ProjectInfo$Modes) & (Mode2%in%PROJECT1$ProjectInfo$Modes))
      {
        POSITION_Mode2<-sum((Mode2==names(PROJECT1$Attributes))*c(1:length(PROJECT1$Attributes)))
        if(sum(colnames(NET1)!=PROJECT1$Attributes[[POSITION_Mode2]]$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The column names of the network are not the same as those for the existing network.\n',
               '##        Check the column names of the network or correctly specify the second value of the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
        #This means that Mode1 (rows) is new
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode1)

        ATTRIBUTEFILE1<-data.frame(rownames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        PROJECT1$Attributes$NEW<-ATTRIBUTEFILE1
        names(PROJECT1$Attributes)[length(PROJECT1$Attributes)]<-Mode1

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode1,paste("Names of the nodes for mode",Mode1)))
      }

      if((!Mode1%in%PROJECT1$ProjectInfo$Modes) & (!Mode2%in%PROJECT1$ProjectInfo$Modes))
      {
        #This means that Mode1 and Mode2 are new
        PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode1,Mode2)

        ATTRIBUTEFILE1<-data.frame(rownames(NET1))
        colnames(ATTRIBUTEFILE1)<-"NodeName"
        ATTRIBUTEFILE2<-data.frame(colnames(NET1))
        colnames(ATTRIBUTEFILE2)<-"NodeName"
        PROJECT1$Attributes$NEW1<-ATTRIBUTEFILE1
        PROJECT1$Attributes$NEW2<-ATTRIBUTEFILE2
        names(PROJECT1$Attributes)[length(PROJECT1$Attributes)-1]<-Mode1
        names(PROJECT1$Attributes)[length(PROJECT1$Attributes)]<-Mode2

        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode1,paste("Names of the nodes for mode",Mode1)))
        PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,c("NodeName",Mode2,paste("Names of the nodes for mode",Mode2)))
      }

      if((Mode1%in%PROJECT1$ProjectInfo$Modes) & (Mode2%in%PROJECT1$ProjectInfo$Modes))
      {
        POSITION_Mode1<-sum((Mode1==names(PROJECT1$Attributes))*c(1:length(PROJECT1$Attributes)))
        POSITION_Mode2<-sum((Mode2==names(PROJECT1$Attributes))*c(1:length(PROJECT1$Attributes)))

        if(sum(rownames(NET1)!=PROJECT1$Attributes[[POSITION_Mode1]]$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The row names of the network are not the same as those for the existing network.\n',
               '##        Check the row names of the network or correctly specify the first value of the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }

        if(sum(colnames(NET1)!=PROJECT1$Attributes[[POSITION_Mode2]]$NodeName)>0)
        {
          stop('\n',
               '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
               '## ERROR: The column names of the network are not the same as those for the existing network.\n',
               '##        Check the column names of the network or correctly specify the second value of the Mode of the new network.\n',
               '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
        }
      }
    }
  }
  PROJECT1[[NLength+1]]<-NET1
  names(PROJECT1)[[NLength+1]]<-NetworkName

  #Add to NetworkInfo
  PROJECT1$NetworkInfo<-rbind(PROJECT1$NetworkInfo,c(NetworkName,Mode1,Mode2,Directed,Loops,Values,Class))
  #Add to NetworksDescription
  PROJECT1$ProjectInfo$NetworksDescription<-rbind(PROJECT1$ProjectInfo$NetworksDescription,c(NetworkName, NetworkDescription))

  cat(' ---------------------------------------------------------------\n\n')
  PROJECT1
}

#' Deletes a network from an existing "dataproject" object
#'
#' Takes a network (specified as "PROJECT"$"NETWORK") as input and deletes it from the overall "dataproject" ("PROJECT") it is part of.
#' This means that the input "NET1" needs to have the form "PROJECT"$"NETWORK" so that the specified network ("NETWORK") can be deleted from the "dataproject" ("PROJECT").
#' - Last updated: 1 May 2022.
#'
#' @param NET1 A network which is part of a "dataproject" (which needs to be of type "PROJECT"$"NETWORK").
#'
#' @return A "dataproject" without the specific network.
#'
#' @export
#'
#' @author Filip Agneessens, \email{filipagneessens2@@gmail.com}
#'
#' @examples
#' ## We can consider an existing network and add the transposed
#' ## Create a network with ordinal categories (0='no relation', 1='acquaintance', 2='friend',
#' ##    3='good friend') as an (adjacency) matrix:
#' Friendship3<-matrix(c(0,1,1, 2,0,3, 0,2,0),3,3, byrow=TRUE)
#' ## Add node names to the matrix:
#' rownames(Friendship3)<-c("a","b","c")
#' colnames(Friendship3)<-rownames(Friendship3)
#'
#' ## Now create the "dataproject", which we call "SchoolClass_Project1"
#' SchoolClass_Project1<-xCreateProject(GeneralDescription="Dataset among 3 students",
#'      NetworkName="Friendship",
#'      NETFILE1=Friendship3,
#'      FileType="Robject",
#'      InFormatType="AdjMat",
#'      NetworkDescription="Strength (0=no relation, 1=acquaintance, 2=friend, 3=good friend)",
#'      Mode=c("Students"),
#'      Directed=TRUE,
#'      Loops=FALSE,
#'      Values="Ordinal",
#'      Class="matrix",
#'      References="No references")
#' ## Check the newly created "dataproject" object
#' SchoolClass_Project1
#'
#' ## Now create the transposed and add:
#' SchoolClass_Project1<-xAddToProject(ProjectName=SchoolClass_Project1,
#'         NetworkName="Friendship_Tr",
#'         NETFILE1=t(SchoolClass_Project1$Friendship),
#'         FileType="Robject",
#'         InFormatType="AdjMat",
#'         NetworkDescription="Transposed of the friendship network",
#'         Mode=c("Students"),
#'         Directed=TRUE,
#'         Loops=FALSE,
#'         Values="Ordinal",
#'         Class="matrix")
#' ## Check the new result:
#' SchoolClass_Project1
#'
#' ## Remove the original friendship network:
#' SchoolClass_Project1<-xRemoveFromProject(SchoolClass_Project1$Friendship)
#' ## Check that network has been removed:
#' SchoolClass_Project1

xRemoveFromProject<-function(NET1)
{
  cat('\n ------ FUNCTION: xRemoveFromProject ---------------------------\n\n')
  # Get name of object NET1 and turn it into string:
  FILENAME1<-deparse(substitute(NET1))
  # Using gregexec the positions where the character in NET1 is "$" are stored in [[1]]
  # The name of the network should be after the last "$", which we can find here:
  MAXDOLLAR<-max(gregexec("\\$",FILENAME1)[[1]])
  # If $ is -1 this means there are no $ in the name, and hence the format is not PROJECT$NETWORK
  # Give error and advice if not at least one $
  if(MAXDOLLAR==-1)
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The input object [', FILENAME1, '] is not a file of format PROJECT$NETWORK', '\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  # Get name of NETWORK
  NETWORKNAME<-substr(FILENAME1,MAXDOLLAR+1,nchar(FILENAME1))
  # Get name of PROJECT (which can be ...$...$... )
  PROJECTNAME<-substr(FILENAME1,1,MAXDOLLAR-1)
  # Provide info of what is attempted:
  cat(' Searching for network [',NETWORKNAME,'] in "dataproject" [',PROJECTNAME,'].\n\n', sep="")

  # Check that "dataproject" exists
  if(class(try(get(PROJECTNAME), silent=TRUE))=="try-error")
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "dataproject" [',PROJECTNAME,'] not found.\n',
         '##        Make sure the name of the "PROJECT" object exists and has the format PROJECT$NETWORK.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }
  # Now retrieve project from memory
  PROJECTOUT<-get(PROJECTNAME)

  # Check if the requested network exist in that project:
  if(!NETWORKNAME %in% names(PROJECTOUT))
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==', '\n',
         '## ERROR: The network [', NETWORKNAME ,'] in [', PROJECTNAME ,'] does not exist.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  # If all fine, remove the network from the existing "dataproject"
  # Get position of network in "dataproject" list object
  POSITION1V<-sum((names(PROJECTOUT)==NETWORKNAME)*c(1:length(PROJECTOUT)))
  # Remove from list object
  PROJECTOUT<-PROJECTOUT[-POSITION1V]
  # Confirm that the network was removed the project
  cat(' Network [',NETWORKNAME,'] removed from "dataproject" [',PROJECTNAME,'].\n\n', sep="")

  # Check if network also in NetworkInfo by checking NetworkName
  if(NETWORKNAME %in% PROJECTOUT$NetworkInfo$NetworkName)
  {
    # Get position of network in NetworkInfo
    POSITION2V<-sum((PROJECTOUT$NetworkInfo$NetworkName==NETWORKNAME)*c(1:dim(PROJECTOUT$NetworkInfo)[1]))
    # Remove row for that network from dataframe
    PROJECTOUT$NetworkInfo<-PROJECTOUT$NetworkInfo[-POSITION2V,]
    # Confirm that the network was removed from PROJECT$NetworkInfo
    cat(' Network [',NETWORKNAME,'] removed from $NetworkInfo in [',PROJECTNAME,'].\n\n', sep="")
  }
  else
  {
    # Give warning that the network did not appear in NetworkInfo
    cat(' *** WARNING: Network [',NETWORKNAME,'] was not found in $NetworkInfo in [',PROJECTNAME,'].\n\n', sep="")
  }
  cat(' ---------------------------------------------------------------\n\n')

  # Check if network also in $ProjectInfo$NetworksDescription by checking NetworkName
  if(NETWORKNAME %in% PROJECTOUT$ProjectInfo$NetworksDescription$NetworkName)
  {
    # Get position of network in NetworkInfo
    POSITION3V<-sum((PROJECTOUT$ProjectInfo$NetworksDescription$NetworkName==NETWORKNAME)*c(1:dim(PROJECTOUT$ProjectInfo$NetworksDescription)[1]))
    # Remove row for that network from dataframe
    PROJECTOUT$ProjectInfo$NetworksDescription<-PROJECTOUT$ProjectInfo$NetworksDescription[-POSITION3V,]
    # Confirm that the network was removed from PROJECT$NetworkInfo
    cat(' Network [',NETWORKNAME,'] removed from $ProjectInfo$NetworksDescription in [',PROJECTNAME,'].\n\n', sep="")
  }
  else
  {
    # Give warning that the network did not appear in NetworkInfo
    cat(' *** WARNING: Network [',NETWORKNAME,'] was not found in $ProjectInfo$NetworksDescription in [',PROJECTNAME,'].\n\n', sep="")
  }
  cat(' ---------------------------------------------------------------\n\n')

  PROJECTOUT
}

#' Adds nodal attributes to an existing "dataproject" object
#'
#' Takes nodal attribute data and adds these to the "$Attributes" part of an existing "dataproject".
#' - Last updated: 1 May 2022.
#'
#' @param ProjectName The name of the "dataproject", to which the attribute needs to be added.
#' @param ATTFILE1 An attribute file that needs to be added. This could be a csv file or an R object (data.frame). The input type is defined in "FileType". If reading in a csv file, the first column should contain the node names (labels), so the function can check that the names correspond with the existing data. If it is an R object, these names should be in the rownames.
#' @param FileType The type of data that is being imported. Default is "csv". Other options are "Robject" when the data is already an object in R, and "csv2". Note that by default a csv file use a point (".") for the decimal point and a comma for the separator in most countries. However, in some Western European countries a comma is used for the decimal point and a semicolon for the separator, which will require "csv2".
#' @param Mode Character (string), indicating the name for the mode of the network.
#' @param AttributesDescription A vector with information about each of the attributes added.
#'
#' @return A "dataproject" object containing the nodal attribute data.
#'
#' @importFrom utils read.csv read.csv2
#'
#' @export
#'
#' @author Filip Agneessens, \email{filipagneessens2@@gmail.com}
#'
#' @examples
#' #' ## We can consider an existing network and add an attribute
#' ## Create a network with ordinal categories (0='no relation', 1='acquaintance', 2='friend',
#' ##    3='good friend') as an (adjacency) matrix:
#' Friendship3<-matrix(c(0,1,1, 2,0,3, 0,2,0),3,3, byrow=TRUE)
#' ## Add node names to the matrix:
#' rownames(Friendship3)<-c("a","b","c")
#' colnames(Friendship3)<-rownames(Friendship3)
#'
#' ## Now create the "dataproject", which we call "SchoolClass_Project1"
#' SchoolClass_Project1<-xCreateProject(GeneralDescription="Dataset among 3 students",
#'      NetworkName="Friendship",
#'      NETFILE1=Friendship3,
#'      FileType="Robject",
#'      InFormatType="AdjMat",
#'      NetworkDescription="Strength (0=no relation, 1=acquaintance, 2=friend, 3=good friend)",
#'      Mode=c("Students"),
#'      Directed=TRUE,
#'      Loops=FALSE,
#'      Values="Ordinal",
#'      Class="matrix",
#'      References="No references")
#' ## Check the newly created "dataproject" object
#' SchoolClass_Project1
#'
#' ## Add nodal attribute data to an existing example "dataproject":
#' Age3<-c(32,65,24)
#' Gender3<-c(1,2,1)
#' X1<-data.frame(Age3,Gender3)
#' rownames(X1)<-c("a","b","c")
#' SchoolClass_Project1<-xAddAttributesToProject(SchoolClass_Project1, ATTFILE1=X1, FileType="Robject",
#'       AttributesDescription=c("Age in years","Gender, where 2=Male and 1=Female"), Mode="Students")
#' ## Check that the attribute data have been added:
#' SchoolClass_Project1

xAddAttributesToProject<-function(ProjectName,
                                  ATTFILE1=NULL,
                                  FileType=NULL,
                                  AttributesDescription=NULL,
                                  Mode=NA)
{
  cat('\n ------ FUNCTION: xAddAttributesToProject ---------------------------\n\n')
  # get project name for object ProjectName
  ProjectName1<-deparse(substitute(ProjectName))

  #check that "dataproject" with name "ProjectName" exists
  if(class(try(get(ProjectName1), silent=TRUE))=="try-error")
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: "ProjectName" not found. Make sure the name of the "dataproject" object exists,\n',
         '##        and does not contain typos.\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }


  # Check if ATTFILE1 is provided
  if(is.null(ATTFILE1)) {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: No ATTFILE1 provided\n',
         '##        Add the file name using [ATTFILE1=...].\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  # Check if FileType is provided
  if(is.null(FileType))
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The FileType should be "csv", "csv2", or "Robject"\n',
         '##        Choose among: [FileType="csv"], [FileType="csv2"], and [FileType="Robject"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  # Check if FileType is a valid option
  if((FileType!="csv") & (FileType!="csv2") & (FileType!="Robject"))
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: The FileType should be "csv", "csv2", or "Robject"\n',
         '##        Choose among: [FileType="csv"], [FileType="csv2"], and [FileType="Robject"]\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  # Check if Mode is provided
  if(is.null(Mode))
  {
    stop('\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Mode should be specified and be a single value of length 1 (e.g.,"Person").\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  # Check if Mode is a single element
  if(length(Mode)!=1)
  {
    stop('\n\n',
         '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
         '## ERROR: Mode should be a vector of length 1 (e.g.,"Person").\n',
         '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
  }

  # Provide info of what is attempted:
  cat(' Getting "dataproject" with name [',ProjectName1,'].\n\n', sep="")

  # now get the "dataproject"
  PROJECT1<-get(ProjectName1)

  # Importing attribute data from a csv file or object
  if(FileType=="csv")
  {
    #check for errors
    ATT1<-read.csv(ATTFILE1, stringsAsFactors=FALSE, row.names=1)
  }

  if(FileType=="csv2")
  {
    #check for errors
    ATT1<-read.csv2(ATTFILE1, stringsAsFactors=FALSE, row.names=1)
  }

  if(FileType=="Robject")
  {
    if(class(ATTFILE1)=="data.frame")
    {
      ATT1<-ATTFILE1
    }
    if(class(ATTFILE1)=="matrix")
    {
      ATT1<-as.data.frame(ATTFILE1)
    }
    else
    {
      if(class(ATTFILE1)!="data.frame")
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: ATTFILE1',ATTFILE1,' should be an object of class data.frame or matrix.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }
    }
  }
  cat(' Attribute file [',deparse(substitute(ATTFILE1)),'] imported\n\n', sep="")
  #-----
  #if Mode==NA
  if (is.na(Mode))
  {
    Mode<-"A"
  }

  #---------------------------------------------------------------------------------------
  #Check names of new file against existing attribute file.
  #if we already have more than one mode in the dataset
  if(length(PROJECT1$ProjectInfo$Modes)>1)
  {
    print(PROJECT1$ProjectInfo$Modes)
    if(Mode %in% PROJECT1$ProjectInfo$Modes)
    {
      POSITION_Mode1<-sum((Mode==names(PROJECT1$Attributes))*c(1:length(PROJECT1$Attributes)))
      if(nrow(ATT1)!=nrow(PROJECT1$Attributes[[POSITION_Mode1]]))
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: The number of rows for the attribute dataset is not the same as that for the exisiting attribute file.\n',
             '##        Check the dataset ATTFILE1.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }
      if(sum(rownames(ATT1)!=PROJECT1$Attributes[[POSITION_Mode1]]$NodeName)>0)
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: The row names of the attribute dataset are not the same as those for the existing attribute file.\n',
             '##        Check the dataset ATTFILE1.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }
      rownames(ATT1)<-NULL
      PROJECT1$Attributes[[POSITION_Mode1]]<-cbind(PROJECT1$Attributes[[POSITION_Mode1]],ATT1)
    }

    if(!(Mode %in% PROJECT1$ProjectInfo$Modes))
    {
      PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode)
      ATTRIBUTEFILE1<-cbind(rownames(ATT1),ATT1)
      colnames(ATTRIBUTEFILE1)[1]<-"NodeName"
      rownames(ATTRIBUTEFILE1)<-NULL
      PROJECT1$Attributes$NEW<-ATTRIBUTEFILE1
      names(PROJECT1$Attributes)[length(PROJECT1$Attributes)]<-Mode

    }
  }

  #One single mode
  if(length(PROJECT1$ProjectInfo$Modes)==1)
  {
    if(Mode==PROJECT1$ProjectInfo$Modes)
    {
      if(nrow(ATT1)!=nrow(PROJECT1$Attributes))
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: The number of rows for the attribute dataset is not the same as that for the exisiting attribute file.\n',
             '##        Check the dataset ATTFILE1.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }
      if(sum(rownames(ATT1)!=PROJECT1$Attributes$NodeName)>0)
      {
        stop('\n',
             '== == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==\n',
             '## ERROR: The row names of the attribute dataset are not the same as those for the existing attribute file.\n',
             '##        Check the dataset ATTFILE1.\n',
             '== == == == Function aborted. Please see error above. == == == == == == == == == == == ==\n')
      }
      rownames(ATT1)<-NULL
      PROJECT1$Attributes<-cbind(PROJECT1$Attributes,ATT1)
    }

    # make new mode-specific attribute
    if(Mode!=PROJECT1$ProjectInfo$Modes)
    {
      PROJECT1$ProjectInfo$Modes<-c(PROJECT1$ProjectInfo$Modes,Mode)
      ATTRIBUTEFILE1<-cbind(rownames(ATT1),ATT1)
      colnames(ATTRIBUTEFILE1)[1]<-"NodeName"
      rownames(ATTRIBUTEFILE1)<-NULL
      PROJECT1$Attributes<-list(PROJECT1$Attributes,ATTRIBUTEFILE1)
      names(PROJECT1$Attributes)<-PROJECT1$ProjectInfo$Modes
    }
  }

  # Adding to $ProjectInfo$AttributesDescription
  ADDATT<-data.frame(colnames(ATT1),Mode,AttributesDescription)
  names(ADDATT)<-c("Variable","Mode","Details")
  PROJECT1$ProjectInfo$AttributesDescription<-rbind(PROJECT1$ProjectInfo$AttributesDescription,ADDATT)

  PROJECT1
}

