#'Demultiplexed data condensing
#'
#'Proccesses demultiplexed fasta files, finds identical reads and deletes them.
#'The number of identical reads is retained in sequence headers.
#'Condensed reads are saved as FASTA file.
#'
#'@param fasta_soubor file demultiplexed data
#'@param prah cutoff for number of unique sequences, sequences with amount less than
#'the cutoff are not saved in the final FASTA file
#'
#'@return None
#'
#'@examples
#'condense("demultiplexed_reads.fasta", prah=10)
#'
#'@export
condense <- function (fasta_soubor, prah=10){

  x<-readDNAStringSet(fasta_soubor)
  if (length(grep("N", x, ignore.case = T))!=0)
  {  
    x <- x[-grep("N", x, ignore.case = T)]
  }
  data<-as.data.frame(BiocGenerics::table(x))
  data<-data[order(data$Freq, decreasing = TRUE), ]
  pocet<-sum(data[,2]>prah)
  kond<-DNAStringSet(as.vector.factor(data[1:pocet,1]))
  pocty<-data[1:pocet,2]
  names(kond)<-paste(pocty,'_sekvenci',sep='')

  nazev<-paste(strsplit(fasta_soubor,'\\.')[[1]][1],'_condensed.fasta',sep='')
  writeXStringSet(kond, nazev)
  return(kond)
}
