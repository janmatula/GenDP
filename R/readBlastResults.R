#'Load results of BLAST database search to R
#'
#'Loads BLAST search results in .txt format to R list structure.
#'
#'@param blastResultsFile .txt file containing results of BLAST database search obtained using the -outfmt 0 command. 
#'@param results number of top results for each query sequence that will be loaded to R
#'
#'@return List containing top BLAST atabase matches for each query sequence.
#'
#'@examples
#'scoreList <- readBlastResults("BLAST_S1_fw_condensed.txt", results = 1000)
#'
#'@export
#'
readBlastResults <- function (blastResultsFile, results = 10)
{
bl <-read.delim(blastResultsFile)
outputPosition <- grep ("Query=", bl[[1]])

scoreMatrix <- matrix(rep(0, 3*results), nrow = results, ncol = 3)
colnames(scoreMatrix)<- c("Name", "BitScore", "Evalue")
scoreList <- replicate(length(outputPosition),scoreMatrix, simplify=FALSE)

for (i in 1:length(outputPosition)){
  names(scoreList)[i]<-substr(bl[[1]][outputPosition[i]],8,nchar(as.character(bl[[1]][outputPosition[i]])))
  
  for (j in 1:results){
    radek <- strsplit(as.character(bl[[1]][j+outputPosition[i]+3]), " ")
    radek<- radek[[1]][radek[[1]] != ""]
    scoreList[[i]][j,1]<- paste(radek[1], radek[2])
    scoreList[[i]][j,2]<- radek[3]
    scoreList[[i]][j,3]<- radek[4]
  }
}
return(scoreList)
}