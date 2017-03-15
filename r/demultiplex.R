#'Raw data demultiplexing
#'
#'User inputs fasta file and MIDs an function outputs a fasta file of demultiplexed reads.
#'Demultiplexed reads are saved as FASTA.
#'
#'@param fasta_soubor file with raw data in FASTA format
#'@param f_primer forward MID
#'@param r_primer reverse MID
#'
#'@return None
#'
#'@examples
#'demultiplex("MHC_reads.fasta", "ACTGACTGGT", "ATGACGTAGT")
#'
#'@export
demultiplex <- function (fasta_soubor, f_primer, r_primer) {

  sek <- readDNAStringSet(fasta_soubor,format='fasta')

  f_primer_revc <- reverseComplement(DNAString(f_primer))
  r_primer_revc <- reverseComplement(DNAString(r_primer))

  r_vyraz_fw <- paste('^',f_primer,'.+',r_primer_revc, '$', sep = "")
  r_vyraz_rev <- paste('^',r_primer,'.+',f_primer_revc, '$', sep = "")


  p_fw <- grep (r_vyraz_fw, sek)
  p_rev <- grep (r_vyraz_rev, sek)

  fw_ready <- sek[p_fw]
  rev_ready <-sek[p_rev]

  fw_ready <- subseq(fw_ready,11,width(fw_ready)-10)
  rev_ready <- subseq(rev_ready,11,width(rev_ready)-10)

  writeXStringSet(fw_ready,'fw_ready.fas')
  writeXStringSet(rev_ready,'rev_ready.fas')

}
