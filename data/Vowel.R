vowel <- read.table("vowel2.dat")
vowel <- as.matrix(vowel)
vowel.mem <- rep(1:11, 48)

Vowel <- list(vowel = vowel, vowel.mem = vowel.mem)

