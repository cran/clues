iris<-read.table("iris.dat", sep=",", header=TRUE)
iris<-as.matrix(iris)
iris.mem<-rep(1:3,rep(50,3))

Iris<-list(iris=iris, iris.mem=iris.mem)

