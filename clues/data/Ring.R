ring<-read.table("broken-ring.dat")
ring<-as.matrix(ring)

# true cluster membership
ring.mem<-scan("broken-ring.mem", quiet=TRUE)

Ring<-list(ring=ring, ring.mem=ring.mem)
