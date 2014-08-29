library("qtl")      # Load the library
data(multitrait)    # Load the data
multitrait          # Show the library

pull.pheno(multitrait, pheno.col = 1)        # Shows the phenotype values of phenotype 1
hist(pull.pheno(multitrait, pheno.col = 1))  # Create a histogram

# Plot a map versus the other ones
plot.map(multitrait, est.map(multitrait))

# Create a histogram using more breaks, add a red line for the cut-off
hist(pull.pheno(multitrait, pheno.col = 1), breaks=100)
abline(v=750, col='red', lty=2, lwd=3)

# Transform out values in a binary phenotype
multitrait$pheno[,1] <- as.numeric(multitrait$pheno[,1] > 750)

npres <- scanone(multitrait, model="np", pheno.col=1)

# Show the differences between binary and non parametric
plot(scanone(multitrait, model="binary", pheno.col=1), npres, col=c("black", "blue"))

data(multitrait) #Reload data, undo the transformation of phenotype 1

# Show where the difference comes from
effectplot(multitrait,mname1="GA1",mname2="GH.117C")

# Permutations
permutationResults <- scanone(multitrait, model ="np", n.perm=1000)

# Map all traits
res <- scanone(multitrait, model="np", pheno.col=1:24)

# Create a heatmap
image(1:117, 1:24, as.matrix(res[,3:26]), breaks =c(0,3,5,10,1000), col=c("white","lightgray","lightblue","black"),xlab="Marker", ylab="Metabolite")
box();grid()
