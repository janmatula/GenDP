#'Transform DNAStringSet to binary representation
#'Transforms DNAStringSet to binary representation
#'@param x sequences as DNAStringSet
#'@return List of sequences in binary represenation.
#'@examples
#'binRep(sequences)
#'@export
binRep<-function(x){
  a<- '00'
  c<- '01'
  g<- '10'
  t<- '11'
  x<-gsub('a', a, x, ignore.case = T )
  x<-gsub('c', c, x, ignore.case = T )
  x<-gsub('g', g, x, ignore.case = T )
  x<-gsub('t', t, x, ignore.case = T )
  return(x)
}
