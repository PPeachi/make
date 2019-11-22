genbank2fasta<-function(file){
  seq<-readLines(file)
  a_n<-paste0(">",unlist(strsplit(seq[1],split = "\\s+"))[2])
  st<-grep("ORIGIN",seq)
  ed<-grep("^//",seq)
  gb_seq<-gsub(" ","",gsub("\\d","",paste(seq[(st+1):(ed-1)],collapse = "")))
  res<-paste0(a_n,"\n",gb_seq)
  return(res)
}
accn<-readLines("accn.txt")
for(i in accn){
  x<-genbank2fasta(paste(i,".gb",sep=""))
  write.table(x,
              file = paste(i,".fas",sep=""),
              row.names = F,
              col.names = F,
              quote = F)
}
