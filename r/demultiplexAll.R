#'Continuous raw data demultiplexing
#'
#'User inputs fasta file and the function then asks him to input MIDs for various reads.
#'Demultiplexed reads are saved as FASTA.
#'
#'@param fasta_soubor file with raw data in FASTA format
#'
#'@return None
#'
#'@examples
#'demultiplex_all("MHC_reads.fasta")
#'
#'@export
demultiplex_all <- function(fasta_soubor){
  sek<-readDNAStringSet(fasta_soubor)

  continue<-'Y'

  while (continue=='Y'){
    fw_primer<- readline('Forward primer:')
    rev_primer<-readline('Reverse primer:')
    nazev<- readline('Nazev readu:')

    fw_primer_revc <- reverseComplement(DNAString(fw_primer))
    rev_primer_revc <- reverseComplement(DNAString(rev_primer))

    r_vyraz_fw <- paste('^',fw_primer,'.+',rev_primer_revc, '$', sep = "")
    r_vyraz_rev <- paste('^',rev_primer,'.+',fw_primer_revc, '$', sep = "")


    fw_ready<- grep(r_vyraz_fw,sek)
    rev_ready<-grep(r_vyraz_rev,sek)

    writeXStringSet(sek[fw_ready], paste(nazev,'_fw_ready.fasta', sep=''))
    writeXStringSet(sek[rev_ready], paste(nazev,'_rev_ready.fasta', sep=''))

    continue<- readline('Pokracovat v zadavani primeru? [Y/N]')
  }
}
