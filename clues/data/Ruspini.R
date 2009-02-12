ruspini<-read.table("ruspini.dat", sep=",", header=TRUE)
ruspini<-as.matrix(ruspini)
ruspini.mem<-c(rep(1,20),rep(2,23),rep(3,17),rep(4,15))

Ruspini<-list(ruspini=ruspini, ruspini.mem=ruspini.mem)

