#'Convert binary to int
#'
#'Tranforms set of binary numbers to integers.
#'
#'@param s set of binary numbers
#'
#'@return None
#'
#'@examples
#'convert(set)
#'
#'@export
convert<-function(s){
  s<-strsplit(s,'')
  f <- function(y) sum(as.numeric(y) * 2^rev((seq_along(as.numeric(y))-1)))
  sapply(s, f)
}
