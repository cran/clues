maronna <- read.table("maronna.dat")
maronna <- as.matrix(maronna)
maronna.mem <- rep(1:4, rep(50,4))

Maronna <- list(maronna = maronna, maronna.mem = maronna.mem)

