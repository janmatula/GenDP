#'Search hash table
#'
#'Searches hash table obtained using the buildHashTable function for occurences of query sequences from FASTA file.
#'
#'@param querySequences FASTA file containing sequences we are searching for
#'@param hashTable output of buildHashTable function
#'@param k number of bases the query sequences are divided into, must be the same as in buildHashTable
#'@param results number of top results of the SSAHA method that will be provie in the outpus list
#'@param tolerance tolerance for insertions and deletions in number of nucleotides
#'
#'@return A list of matrices containg best matches and score for each query sequence.
#'
#'@examples
#'searchHashTable("query.fasta", hashTable, k=10, results = 10, tolerance = 5)
#'
#'@export
searchHashTable <- function(querySequences, hashTable, k = 10, results = 1, tolerance=1)
{
L <- hashTable[[3]]
A <- hashTable[[2]]
query <- readDNAStringSet(querySequences)
query <- binRep(query)
query <- lapply(query, function(query) sapply(seq(from = 1, to = nchar(query) - k * 2 + 1, by = 2), function(i) convert(substr(query, i, i + k * 2 - 1))))

sequence_identification_matrix <- matrix(rep(0), 3 * results, nrow = results, ncol=3)
colnames(sequence_identification_matrix) <- c("Db number","Name", "Score")
sequence_identification <- replicate(length(query),sequence_identification_matrix, simplify = FALSE)
names(sequence_identification)<-paste0(names(query))

for (ind in 1:length(query)) {
  p <- sapply(query[[ind]], function(query) match(query, A))
  if (all(is.na(p) == TRUE)) {
    next
  }
  M <- L[p]
  t <- 0:length(M)
  M <- lapply(1:length(M), function(i) cbind(M[[i]], M[[i]][, 2] - t[i]))
  M <- M[sapply(M, function(x) dim(x)[1]) > 0]
  M <- do.call(rbind, M)
  M[, c(1, 2, 3)] <- M[, c(1, 3, 2)]
  M <- matrix(M[order(M[, 1], M[, 2]), ], ncol = 3)
  
  
  u <- unique(M[, 1])
  score <- rep(0, length(u))
  
  for (j in 1:length(u)) {
    o <- which(M[, 1] == u[j])
    data <- as.data.frame(table(M[o, 2]))
    data[,2]<-data[,2]**2
    if(length(data[,1])>1){
      for (k in 1:(length(data[,1])-1)){
        if (as.numeric(data[k,1])-as.numeric(data[k+1,1])<=tolerance){
          data[k+1,2]<-data[k,2]+data[k+1,2]
        }
      }
    }
    score[j] <- max(data[, 2])
  }
  for (l in 1:results){
    sequence_identification[[ind]][l,1] <- u[which.max(score)]
    sequence_identification[[ind]][l,2] <- hashTable[[1]][u[which.max(score)]]
    sequence_identification[[ind]][l,3] <- max(score)
    score[which.max(score)]<-1
  }
}
return(sequence_identification)
}