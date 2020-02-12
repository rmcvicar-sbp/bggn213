
#Can you improve this analysis code?
#original code from lecture6 worksheet

s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")

#let's reduced it to 1 for simplicity sake 
s1 <- read.pdb("4AKE") 
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")

#Me trying to break this beast down
#step1
s_pdbID <- function(x) {
  read.pdb(file = x)
}
#step2
s_pdbID_chainA <- function(x) {
  trim.pdb((s_pdbID(x)), chain="A", elety="CA")
}
#step3
s_pdbID_chainA.b <- s_pdbID_chainA$atom$b

#step 4, this would be awesome if it worked
plot_HW6 <- function(x) {
plotb3(s_pdbID_chainA.b, sse=s_pdbID_chainA, typ="l", ylab="Bfactor")
}

####not too many issues, let's mess with it
#input for each component
# save each element that we want in our equation
s.pdb <- read.pdb((file[i]))
s.pdb.chainA <- trim.pdb(s.pdb, chain="A", elety="CA")
s.pdb.chainA.b <- s.pdb.chainA$atom$b

#tie it all together....
plot.HW6 <- function(x) {
  plotb3(s.pdb.chainA.b, sse=s.pdb.chainA, typ="l", ylab="Bfactor")
}
###not working, 
#function (x, start, stop) 
#if (!is.character(x)) 
    #x <- as.character(x)
  #Internal(substr(x, as.integer(start), as.integer(stop)))

plot.HW6 <- function(file, chain="A", elety="CA", atom$b) 
  for (i in 1:length(file)) {
    s.pdb <- read.pdb((file[i]))
    s.pdb.chainA <- trim.pdb(s.pdb, chain="A", elety="CA")
    s.pdb.chainA.b <- s.pdb.chainA$atom$b
  }
  
