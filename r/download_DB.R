#'Database download
#'
#'Downloads database from EBI ftp server.
#'
#'@param URL URL to database on EBI
#'@param type type of sequences (nuc for nucletide, prot for aminoacid)
#'@param directory directory to which the database is saved
#'@param mergeIntoOneFasta TRUE: merge all FASTA files to a single FASTA file
#'@param name name of the merged FASTA file
#'
#'@return None
#'
#'@examples
#'download_DB('ftp://ftp.ebi.ac.uk/pub/databases/ipd/mhc/nhp/')
#'
#'@export
download_DB <- function(URL, type='nuc', directory=getwd(), mergeIntoOneFasta=FALSE, name='database'){

  download.file(URL, 'html.csv')
  tabulka <- read.table("html.csv" , fill = TRUE ,row.names=NULL)
  file.remove('html.csv')
  tabulka<-tabulka[5]
  tabulka<-unlist(tabulka)
  r_vyraz<- paste('.*_', type, '.fasta"', sep='')

  linky <- sub(r_vyraz,'\\1',tabulka[grepl(r_vyraz, tabulka)])
  linky <- gsub('</A>', '',linky)
  linky <- gsub('>', URL, linky)

  stara_wd<-getwd()
  setwd(directory)

  sek<-DNAStringSet()
  if (mergeIntoOneFasta==FALSE){
    for (i in 1:length(linky)){
      download.file(linky[i], substring(linky[i], width(URL)+1))
    }
  }
  else {
    for (i in 1:length(linky)){
      download.file(linky[i], substring(linky[i], width(URL)+1))
      sek<-DNAStringSet(c(sek,readDNAStringSet(substring(linky[i], width(URL)+1))))
      file.remove(substring(linky[i], width(URL)+1))
    }
  }
  writeXStringSet(sek, paste(name,'.fasta'))

  setwd(stara_wd)
}
