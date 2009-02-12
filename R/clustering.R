# clustering
clustering<-function(y, disMethod="Euclidean")
{
  disMethod=match.arg(disMethod, c("Euclidean", "1-corr"))

  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  n<-as.integer(nrow(y))
  p<-as.integer(ncol(y))

  n1<-n-1

  if(disMethod == "Euclidean") {
    disMethod2=1
  } else { 
    disMethod2=2
  } 

  #input
  storage.mode(y)<-"double"
  storage.mode(n)<-"integer"
  storage.mode(n1)<-"integer"
  storage.mode(p)<-"integer"
  storage.mode(disMethod2)<-"integer"

  #output
  point<-rep(0, n)
  storage.mode(point)<-"integer"
  db<-rep(0.0, n)
  storage.mode(db)<-"double"
  omin<-0.0
  storage.mode(omin)<-"double"
  nClusters<-0
  storage.mode(nClusters)<-"integer"
  mem<-rep(0.0, n)
  storage.mode(mem)<-"integer"
  size<-rep(0.0, n)
  storage.mode(size)<-"integer"
  
  res<-.Fortran("clustering", y, n, n1, p, disMethod2, point=point, 
    db=db, omin=omin, nClusters=nClusters, mem=mem, size=size, 
    PACKAGE="clues") 

  g<-res$nClusters
  size<-res$size[1:g]
  db<-res$db[1:n1]
  resList<-list(mem=res$mem, size=size, g=g, 
    db=db, point=res$point, omin=res$omin)

  return(resList)
}
