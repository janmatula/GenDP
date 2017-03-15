#'Build hash dable
#'
#'Builds hash table from FASTA file containing database.
#'
#'@param fasta_soubor FASTA file containing database
#'@param k number of bases the database is divided into
#'
#'@return List containing list of pointers, list of indexes and list of sequence headers
#'
#'@examples
#'buildHashTable("database.fasta", k=10)
#'
#'@export
buildHashTable<- function (fasta_soubor, k=10){
  sek<-readDNAStringSet(fasta_soubor)
  ids<-names(sek)
  sek<-gsub('N', 'A', sek)
  sek<-binRep(sek)
  p<-list()

  p<-lapply(sek, function(sek) sapply(seq(from=1, to=nchar(sek)-k*2+1, by=2*k), function(i) substr(sek, i, i+k*2-1)))
  p<-lapply(p, function(p) convert(p))

  data<-as.data.frame(table(unlist(p)))
  A<-as.numeric(as.vector.factor(data[,1]))

  L<-lapply(1:length(A), function(x) matrix(0, nrow=data[,2][x], ncol=2) )
  pocitadlo<-rep(1,length(A))

  for (i in 1:length(p)){
    pp<-1
    for (j in 1:length(p[[i]])){
      poz<-which(A==p[[i]][j])
      L[[poz]][pocitadlo[poz],1]<-i
      L[[poz]][pocitadlo[poz],2]<-pp
      pocitadlo[poz]<- pocitadlo[poz]+1
      pp<-pp+k
    }

  }
  return(list(ids, A, L))
}
