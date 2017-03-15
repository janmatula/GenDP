#'Search hash table
#'
#'Searches hash table form buildHashTable for query sequences from FASTA file.
#'
#'@param querySequences FASTA file containing sequences we are searching for
#'@param hashTable output of buildHashTable function
#'@param k number of bases the query sequences are divided into, must be the same as in buildHashTable
#'
#'@return A matrix containg best match and score for each sequence
#'
#'@examples
#'searchHashTable("query.fasta", hashTable k=10)
#'
#'@export
searchHashTable<-function(querySequences, hashTable, k=10){
  L<-hashTable[[3]]
  A<-hashTable[[2]]
  query<-readDNAStringSet(querySequences)
  query<-binRep(query)

  query<-lapply(query, function(query) sapply(seq(from=1, to=nchar(query)-k*2+1, by=2), function(i) convert(substr(query, i, i+k*2-1))))

  sequence_identification<-matrix(rep(0, 2*length(query)), nrow=2)
  rownames(sequence_identification)<-c('Seq Number', 'Score')

  for (ind in 1:length(query)){
    p<-sapply(query[[ind]], function(query) match(query,A))
    if (all(is.na(p)==TRUE)) {next}
    M<-L[p]

    t<-0:length(M)
    M<-lapply(1:length(M), function(i) cbind(M[[i]], M[[i]][,2]-t[i]))
    M<-M[sapply(M, function(x) dim(x)[1]) > 0]
    M<-do.call(rbind,M)
    M[,c(1,2,3)]<- M[,c(1,3,2)]

    M<-matrix(M[order(M[,1],M[,2]),], ncol=3)

    u<- unique(M[,1])

    score<-rep(0, length(u))
    for (j in 1:length(u)){
      o<-which(M[,1]==u[j])
      data<-as.data.frame(table(M[o,2]))
      score[j]<-max(data[,2])
    }

    sequence_identification[1,ind]<-which.max(score)
    sequence_identification[2,ind]<-max(score)
  }
  return(sequence_identification)
}
