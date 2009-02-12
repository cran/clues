# using silhouette index as the measure of the strength of the clusters
clues.sil<-function(y, n0=5, alpha=0.05, eps=1.0e-4, itmax=500,  
  K2.vec=n0, s2=-3, disMethod="Euclidean", 
  plotFlag=FALSE, plot.dim=c(1,2), quiet=FALSE)
{
  disMethod=match.arg(disMethod, c("Euclidean", "1-corr"))

  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  # first pass
  second<-FALSE
  res<-ChooseK.sil(y, y2=y, n0, alpha, eps, itmax, second, K2.vec, s2, 
    disMethod, plotFlag, plot.dim, quiet)

#  if(!quiet)
#  { cat("K.vec>>\n")
#    print(res$K.vec)
#    cat("g.vec>>\n")
#    print(res$g.vec)
#  }
  # check if we need second pass
  # if length(g.vec) = 1 or g.vec[pos-1]=g.vec[pos]+1=g.vec[pos+1]-1,
  # then do NOT need second pass
  len<-length(res$g.vec)
  if(len>1) # may need second pass
  { pos<-which(res$g.vec == res$g)
#    if(!quiet)
#    {
#      cat("final.g=", res$g,"\n")
#      cat("pos>>\n")
#      print(pos)
#    }
    len2<-length(pos)
    len3<-length(res$g.vec)
    pos1<-pos[1]
    pos2<-pos[len2]
    g<-res$g
    s2<-res$avg.s
    K12<-res$K.vec[pos1]
    K21<-res$K.vec[pos2]
    if(pos1>1)
    { g1<-res$g.vec[pos1-1]
      K11<-res$K.vec[pos1-1]
#      if(!quiet)
#      {
#        cat("pos1=",pos1, " g1=", g1, "\n")
#      }
    }
    if(pos2<len3)
    { K22<-res$K.vec[pos2+1]
      g2<-res$g.vec[pos2+1]
#      if(!quiet)
#      { cat("pos2=", pos2, " g2=", g2, "\n") }
    }
    if(pos1>1 && g1-g>1) # need second pass
    { K2.vec<-c(K11, K12)
#      if(!quiet)
#      { cat("K2.vec>>", K2.vec,"\n")
#      }
      # second pass
      second<-TRUE
      res2<-ChooseK.sil(y, y2=res$y.old1, n0, alpha, eps, itmax, second, 
        K2.vec, s2, disMethod, plotFlag, plot.dim, quiet)
      if(res2$myupdate == TRUE) # update the final partition
      { s2<-res2$avg.s 
        res<-res2
      }
    }
    if(pos2<len3 && g-g2>1) # need second pass
    { K2.vec<-c(K21, K22)
#      if(!quiet)
#      { cat("K2.vec>>", K2.vec,"\n")
#      }
      # second pass
      second<-TRUE
      res3<-ChooseK.sil(y, y2=res$y.old2, n0, alpha, eps, itmax, second, 
        K2.vec, s2, disMethod, plotFlag, plot.dim, quiet)
      if(res3$myupdate == TRUE) { res<-res3 }
    }
  }
  return(res)
}

# first pass
# if second=1, then it is the first pass
# if second=2, then it is the second pass
# y --- original data
# y2 --- data used to do shrinking and clustering
ChooseK.sil<-function(y, y2, n0=5, alpha=0.05, eps=1.0e-4, itmax=500, 
  second=F, K2.vec=n0, s2=-3, disMethod="Euclidean", 
  plotFlag=FALSE, plot.dim=c(1,2), quiet=FALSE)
{
  disMethod=match.arg(disMethod, c("Euclidean", "1-corr"))

  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  # step 1
  y<-as.matrix(y)
  y2<-as.matrix(y2)
  dat<-y
  dat2<-y2
  nObs<-nrow(dat)
  nVars<-ncol(dat)
  nClusters0<-n0

  if(disMethod == "Euclidean") {
    disMethod2=1
  } else { 
    disMethod2=2
  } 

  # input
  nObs1<-nObs-1
  storage.mode(dat)<-"double"
  storage.mode(dat2)<-"double"
  storage.mode(nObs)<-"integer"
  storage.mode(nObs1)<-"integer"
  storage.mode(nVars)<-"integer"
  storage.mode(nClusters0)<-"integer"
  storage.mode(alpha)<-"double"
  storage.mode(eps)<-"double"
  ITMAX<-itmax
  storage.mode(ITMAX)<-"integer"
  storage.mode(second)<-"logical"
  nNeiVec2<-K2.vec
  storage.mode(nNeiVec2)<-"integer"
  storage.mode(s2)<-"double"
  storage.mode(disMethod2)<-"integer"
  storage.mode(quiet)<-"logical"
  
  # output
  avgsFinal<-0
  storage.mode(avgsFinal)<-"double"
  sFinal<-rep(0, nObs)
  storage.mode(sFinal)<-"double"
  memFinal<-rep(0, nObs)
  storage.mode(memFinal)<-"integer"
  nClustersFinal<-0
  storage.mode(nClustersFinal)<-"integer"
  clustSizeFinal<-rep(0, nObs)
  storage.mode(clustSizeFinal)<-"integer"
  nNeiFinal<-0
  storage.mode(nNeiFinal)<-"integer"
  nClustVec<-rep(0, nObs)
  storage.mode(nClustVec)<-"integer"
  nNeiVec<-rep(0, nObs)
  storage.mode(nNeiVec)<-"integer"
  myt<-0
  storage.mode(myt)<-"integer"
  datold1<-matrix(0, nrow=nrow(dat), ncol=ncol(dat))
  datold2<-matrix(0, nrow=nrow(dat), ncol=ncol(dat))
  storage.mode(datold1)<-"double"
  storage.mode(datold2)<-"double"
  myupdate<-FALSE
  storage.mode(myupdate)<-"logical"

  res<-.Fortran("chooseksil", dat, dat2, nObs, nObs1, nVars, nClusters0, 
     alpha, eps,
     ITMAX, second, nNeiVec2, s2,  disMethod2, quiet,
     avgsFinal=avgsFinal, 
     sFinal=sFinal,
     memFinal=memFinal, 
     nClustersFinal=nClustersFinal, 
     clustSizeFinal=clustSizeFinal, 
     nNeiFinal=nNeiFinal, 
     nClustVec=nClustVec, 
     nNeiVec=nNeiVec, 
     myt=myt, 
     datold1=datold1, 
     datold2=datold2, 
     myupdate=myupdate,
     PACKAGE="clues")
     resList<-list(K=res$nNeiFinal, 
     size=res$clustSizeFinal[1:res$nClustersFinal],
     mem=res$memFinal, 
     g=res$nClustersFinal, 
     avg.s=res$avgsFinal,
     s=res$sFinal, 
     K.vec=res$nNeiVec[1:res$myt], 
     g.vec=res$nClustVec[1:res$myt],
     myupdate=res$myupdate, 
     y.old1=res$datold1,
     y.old2=res$datold2)

  if(!quiet)
  {
#    cat("final g=", resList$g, " K=", resList$K, "\n")
    if(plotFlag)
    { plotClusters(y, resList$mem, plot.dim) }
  }

  return(resList)
}

Get.Silhouette<-function(y, mem, disMethod="Euclidean")
{
  disMethod=match.arg(disMethod, c("Euclidean", "1-corr"))

  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  n<-nrow(y)
  p<-ncol(y)

  g<-length(unique(mem))
  size<-tapply(mem, mem, length)

  if(disMethod == "Euclidean") {
    disMethod2=1
  } else { #if (disMethod == "corr") {
    disMethod2=2
  } 

  #input
  storage.mode(y)<-"double"
  storage.mode(n)<-"integer"
  storage.mode(p)<-"integer"
  storage.mode(g)<-"integer"
  storage.mode(size)<-"integer"
  storage.mode(mem)<-"integer"
  storage.mode(disMethod2)<-"integer"

  # output
  avgs<-0
  storage.mode(avgs)<-"double"
  memNei<-rep(0, n)
  storage.mode(memNei)<-"integer"
  sIndex<-rep(0, n)
  storage.mode(sIndex)<-"double"

  res<-.Fortran("silhouette", y, n, p, mem, g, size, disMethod2,
    memNei=memNei, sIndex=sIndex, avgs=avgs, PACKAGE="clues") 

  avgs<-res$avgs
  s<-res$sIndex
  neighbor<-res$memNei

  res<-list(avg.s=avgs, s=s, neighbor=neighbor)
}


