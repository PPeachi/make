download_genbank <- function(accession) {
  for (acc in accession) {
    URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", paste(acc, collapse = ","), "&rettype=gb&retmode=text",sep = "")
    #utils::download.file(url = URL, destfile = paste0(acc, ".gb"), quiet = TRUE)
    cmd = paste('curl', paste0("\'", URL, "\'"), '-o', paste0(acc, ".gb"))
    print(cmd)
    system(cmd)
  }
}
acc1<-"AJ534526"
download_genbank(acc1)
