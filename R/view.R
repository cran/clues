# plot trajectory for each cluster
plotCurves<-function(y, mem,
                     myxlab="variable",
                     myylab="observation"
                     )
{
  dat<-y

  memGenes<-mem
  nClusters<-length(unique(mem))
  size<-tapply(mem, mem, length)

  nr<-nrow(dat) # number of observations
  nc<-ncol(dat) # number of variables
 
  ylim<-range(as.vector(dat))
  set1<-1:nc

  if(nClusters<3) # 1 or 2
  { 
    nPanels<-nClusters  
    par(mfrow=c(1, nPanels))
  } else if (nClusters <5) { # 3 or 4
    nPanels<-2 
    par(mfrow=c(2, 2))
  } else {
    nPanels<-ceiling(nClusters/3)
    par(mfrow=c(nPanels, 3))
  }

  for(i in 1:nClusters)
  { dati<-dat[memGenes == i,,drop=FALSE] 
    plot(set1, dati[1,], xlab=myxlab, ylab=myylab, 
      xlim=c(0,nc+1), ylim=ylim, 
      type="l", lty=i, col=i, lwd=3, axes=FALSE) 
    if(size[i]>1)
    { for(j in 2:size[i])
      { lines(set1, dati[j,], lty=i, col=i, lwd=3) }
    }
    axis(1, at=1:nc, labels=1:nc)
    axis(2)
    box()
    title(main="", sub=paste("cluster", i))
  }
  par(mfrow=c(1, 1))

}

# plot average trajectories for each cluster
plotAvgCurves<-function(y, mem,
                     myxlab="variable",
                     myylab="average observation"
                     )
{
  dat<-y

  memGenes<-mem
  nClusters<-length(unique(mem))

  nr<-nrow(dat) # number of observations
  nc<-ncol(dat) # number of variables
 
  ylim<-range(as.vector(dat))
  set1<-1:nc

  if(nClusters<3) # 1 or 2
  { 
    nPanels<-nClusters  
    par(mfrow=c(1, nPanels))
  } else if (nClusters<5) { # 3 or 4
    nPanels<-2 
    par(mfrow=c(2, 2))
  } else {
    nPanels<-ceiling(nClusters/3)
    par(mfrow=c(nPanels, 3))
  }

  for(i in 1:nClusters)
  { dati<-dat[memGenes == i,,drop=FALSE] 
    meancol<-apply(dati, 2, mean)
    plot(set1, meancol, xlab=myxlab, ylab=myylab, 
    xlim=c(0,nc+1), ylim=ylim, 
    type="l", lty=i, col=i, lwd=3, axes=FALSE)
    lines(set1, meancol, lty=i, col=i, lwd=3) 
    axis(1, at=1:nc, labels=1:nc)
    axis(2)
    box()
    title(main="", sub=paste("cluster", i))
  }
  par(mfrow=c(1, 1))
}



plotClusters<-function(y, mem, plot.dim=c(1,2))
{
  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  if(ncol(y)>1)
  {
    umem<-unique(mem)
    
    g<-length(unique(mem))
    t1<-plot.dim[1]
    t2<-plot.dim[2]
    y1<-as.vector(y[, t1])
    y2<-as.vector(y[, t2])
    myxlim<-range(y1)
    myylim<-range(y2)
    mylim<-c(min(myxlim[1], myylim[1]), max(myxlim[2], myylim[2]))
  
    plot(y[mem == umem[1], plot.dim], xlim=mylim, ylim=mylim, 
      xlab=paste("dim",t1),ylab=paste("dim", t2), col=umem[1], pch=umem[1])
    if(g>1)
    { for(i in 2:g)
      { if(sum(mem == umem[i]) == 1)
        { points(y[mem == umem[i],plot.dim[1]], y[mem == umem[i], plot.dim[2]], col=umem[i],pch=umem[i]) }
        else { points(y[mem == umem[i],plot.dim], col=umem[i],pch=umem[i]) }
      }
    }
  } else {
    print("Warning: y is a vector, not a matrix. No plot is outputed")
  }
}


