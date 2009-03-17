# shrinking clusters
shrinking<-function(y, K, disMethod="Euclidean", eps=1.0e-4, itmax=20)
{
  disMethod=match.arg(disMethod, c("Euclidean", "1-corr"))

  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  y<-as.matrix(y)
  n<-nrow(y)
  p<-ncol(y)

  if(disMethod == "Euclidean") {
    disMethod2=1
  } else { #if (disMethod == "corr") {
    disMethod2=2
  } 

  K2<-K+1

  #input
  storage.mode(y)<-"double"
  storage.mode(n)<-"integer"
  storage.mode(p)<-"integer"
  storage.mode(K)<-"integer"
  storage.mode(K2)<-"integer"
  storage.mode(eps)<-"double"
  storage.mode(itmax)<-"integer"
  storage.mode(disMethod2)<-"integer"

  #output
  ynew<-matrix(0, nrow=n, ncol=p)
  storage.mode(ynew)<-"double"
 
  res<-.Fortran("sharpen", y, n, p, K2, K, itmax, eps,
    disMethod2, ynew=ynew, PACKAGE="clues") 

  y.new<-res$ynew
  invisible(y.new)
}


