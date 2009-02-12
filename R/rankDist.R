rankCorr = function(dat, method) {
 cor(t(dat), method=method)
}

rankDist<-function(dat, method="spearman")
{
  if(!is.matrix(dat))
  { dat<-matrix(dat, ncol=1) }
  mat<-rankCorr(dat, method)
  dMat<-as.dist((1-mat)/2)
  return(dMat)
}

  
getDisti<-function(i, dat, method="euclidean")
{

  if(!is.matrix(dat))
  { dat<-matrix(dat, ncol=1) }
  i<-as.integer(i)
  n<-as.integer(nrow(dat))
  p<-as.integer(ncol(dat))
  dat2<-as.double(as.vector(t(dat)))
  tt<-.C("getDisti", i, dat2, n, p, di=as.double(rep(0, n)), 
             PACKAGE="clues")
  di<-tt$di
  return(di)
}

getDistij<-function(i, j, dat, method="euclidean")
{

  if(!is.matrix(dat))
  { dat<-matrix(dat, ncol=1) }
  i<-as.integer(i)
  j<-as.integer(j)
  p<-as.integer(ncol(dat))
  dat2<-as.double(as.vector(t(dat)))
  tt<-.C("getDistij", i, j, dat2, p, dij=as.double(0), 
             PACKAGE="clues")
  dij<-tt$dij
  return(dij)
}

getDistSets<-function(set1, set2, dat, method="euclidean")
{

  if(!is.matrix(dat))
  { dat<-matrix(dat, ncol=1) }
  set1<-as.integer(set1)
  set2<-as.integer(set2)
  n1<-as.integer(length(set1))
  n2<-as.integer(length(set2))
  p<-as.integer(ncol(dat))
  dat2<-as.double(as.vector(t(dat)))
  tt<-.C("getDistSets", set1, n1, set2, n2, dat2, p, dMat=as.double(rep(0, n1*n2)), 
             PACKAGE="clues")
  dMat<-matrix(tt$dMat, nrow=n1, ncol=n2, byrow=T)
  return(dMat)
}

