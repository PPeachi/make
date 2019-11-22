read_fas<-function(file){
  line<-readLines(file)
  i<-grep(">",line)
  id<-sub(">","",line[i])
  start<-(i+1)
  end<-c(i[-1]-1,length(line))
  seq<-sapply(seq_along(start), function(i) paste0(line[start[i]:end[i]],collapse = ""))
  df<-data.frame(id=id,seq=seq)
  df$id<-as.character(df$id)
  df$seq<-as.character(df$seq)
  return (df)
}
x1<-read_fas("AJ534526.fas")
x2<-read_fas("AJ534527.fas")
x<-x1$seq
y<-x2$seq
global_align<-function(x,y,match,mismatch,gap){
  m<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  n<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  l<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  for(j in 2:ncol(m)){
    m[1,j]<-(m[1,j-1]+gap)
    n[1,j]<-'←'
    l[1,j]<-'left'
  }
  for(i in 2:nrow(m)){
    m[i,1]<-(m[i-1,1]+gap)
    n[i,1]<-'↑'
    l[i,1]<-'up'
  }
  for(i in 2:nrow(m)){
    for(j in 2:ncol(m)){
      xl<-unlist(strsplit(x,split = ""))
      yl<-unlist(strsplit(y,split = ""))
      if(xl[i-1]==yl[j-1]){s1<-m[i-1,j-1]+match} else {s1<-m[i-1,j-1]+mismatch}
      s2<-(m[i-1,j]+gap)
      s3<-(m[i,j-1]+gap)
      m[i,j]<-max(s1,s2,s3)
      if(s1==m[i,j]){n[i,j]<-'↖';l[i,j]<-'diag'} 
      if(s2==m[i,j]){n[i,j]<-'↑';l[i,j]<-'up'}
      if(s3==m[i,j]){n[i,j]<-'←';l[i,j]<-'left'}
      if(m[i,j]==s1&&m[i,j]==s2){n[i,j]<-paste0('↖','↑');l[i,j]<-'diag_up'}
      if(m[i,j]==s1&&m[i,j]==s3){n[i,j]<-paste0('↖','←');l[i,j]<-'diag_left'}
      if(m[i,j]==s2&&m[i,j]==s3){n[i,j]<-paste0('↑','←');l[i,j]<-'up_left'}
      if(m[i,j]==s1&&m[i,j]==s2&&m[i,j]==s3){n[i,j]<-paste0('↖','↑','←');l[i,j]<-'diag_up_left'}
    }
  }
  cat("Sequence x:",x,"\n")
  cat("Sequence y:",y,"\n")
  cat("Scoring system:",match,"for match,",mismatch,"for mismatch,",gap,"for gap","\n","\n")
  cat("Dynamic programming matrix 1:","\n")
  print(m)
  cat("Dynamic programming matrix 2:","\n")
  print(n)
  cat("Dynamic programming matrix 3:","\n")
  print(l)
  writeLines(paste0("\n","Alignment:"))
  lx<-unlist(strsplit(x,split = ""))
  ly<-unlist(strsplit(y,split = ""))
  L2=length(lx)
  L1=length(ly)
  D<-c('left','up','diag')
  if(l[L2+1,L1+1]==D[1]){ 
    RX2<-c(ly[L1]);RX1<-c("-");i=L2+1;j=L1 
  } else {
    if(l[L2+1,L1+1]==D[2]){
      RX2<-c("-");RX1<-c(lx[L2]);i=L2;j=L1+1
    } else { 
      RX2<-c(ly[L1]);RX1<-c(lx[L2]);i=L2;j=L1 
    } 
  } 
  while((i>1)&&(j>1)){ 
    # browser() 
    if(l[i,j]==D[1]){ 
      RX2<-c(ly[j-1],RX2);RX1<-c("-",RX1);j=j-1 
    } 
    else if(l[i,j]==D[2]){ 
      RX2<-c("-",RX2);RX1<-c(lx[i-1],RX1);i=i-1 
    } 
    else {RX2<-c(ly[j-1],RX2);RX1<-c(lx[i-1],RX1);j=j-1;i=i-1} 
  }
  RX3<-c()
  for(r in 1:length(RX1)){
    if(RX1[r]==RX2[r]){RX3<-c(RX3,"|")}
    else {RX3<-c(RX3," ")}
  }
  #hamming_distance<-(max(nchar(x),nchar(y))-length(grep("[|]",RX3)))
  if(length(lx)>length(ly)){lmer<-ly;s<-lx} else {lmer<-lx;s<-ly}
  hamming_distance<-numeric(length(s)-length(lmer)+1)
  for(i in 1:(length(s)-length(lmer)+1)){
    ss<-s[i:(i+length(lmer)-1)]
    for(j in 1:length(lmer)){
      if(lmer[j]!=ss[j]){hamming_distance[i]<-(hamming_distance[i]+1)}
    }
  }
  hamming_distance<-min(hamming_distance)
  cat(" x: ",RX1,"\n","   ",RX3,"\n","y: ",RX2,"\n")
  writeLines(paste0("\n","#1 score: ",m[L2+1,L1+1],"\n","#2 hamming-distance: ",hamming_distance))
}
z<-global_align(paste(unlist(strsplit(x,split = ""))[1:10], collapse = ""),paste(unlist(strsplit(y,split = ""))[1:10], collapse = ""),1,-1,-3)
