curve <- read.table("curve.dat")
curve <- as.matrix(curve)

# true cluster membership
curve.mem <- scan("curve.mem", quiet = TRUE)

Curve <- list(curve = curve, curve.mem = curve.mem)

