#'Compare results of SSAHA database search to results of the same BLAST atabase search
#'
#'Compares top results of SSAHA database search to top results of BLAST database search.
#'Outputs the DSC parameter which can be used to test the validity of SSAHA search results as compared to BLAST.
#'
#'@param SSAHAresults  list of top 10 results of obtained by using the searchHashTableFunction.
#'@param blastResultsFile .txt file containing results of BLAST database search obtained using the -outfmt 0 command.
#'@param results number of top results of the SSAHA metod for each query sequence that will be compared to BLAST search results
#'
#'@return Vectorcontaining the DSC percentage parameter for each query sequence.
#'
#'@examples
#'DSC <- compareToBlast(s11_rev_SSAHAresults, "Blast_s11_rev.txt", results =10)
#'
#'@export
#'
compareToBlast <- function (SSAHAresults, blastResultsFile, results = 10)
{
  scoreList <- readBlastResults(blastResultsFile, results = 1000)
  
  percentages <- rep(0, length(x))
  for (i in 1:length(x)){
    v<-rep(0, results)
    ssahaScores<-unique(x[[i]][,3])
    blastScores<-unique(scoreList[[i]][,2])
    scoringVector<-rep(length(blastScores),results)
    for (j in 1:results){
      if (length(which(blastScores == scoreList[[i]][,2][match(as.character(x[[i]][j,2]),scoreList[[i]][,1])]))>0){
        v[j] <- which(blastScores == scoreList[[i]][,2][match(as.character(x[[i]][j,2]),scoreList[[i]][,1])])
      }
    }
    u <- sort(v)
    scoringVector[u==0]<-0
    scoringVector[v==0]<-0
    u[u>0] <- u[u>0]-(min(u[u>0])-1)
    
    w <- unique(u[u>0])
    if(length(w)>1){
      for (k in 1:(length(w)-1)){
        differ <- w[k+1]-w[k]
        if (differ != 1){
          scoringVector<-scoringVector * ((length(blastScores)-differ)/length(blastScores))
        }
      }
    }
    scoringVector <- scoringVector-abs(v-u)
    scoringVector[scoringVector<0]<-0
    percentages[i] <- sum(scoringVector)/(results*length(blastScores))*100
  }
  return(percentages)
}