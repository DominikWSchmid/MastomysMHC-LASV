
########## Code for MHC Supertyping in A. jamaicensis, July 2020
## Ulm University, Institute of Evolutionary Ecology and Conservation Genomics
## Ramona Fleischer



## useful references

# https://academic.oup.com/bioinformatics/article/24/11/1403/191127
# https://grunwaldlab.github.io/Population_Genetics_in_R/clustering_plot.html
# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf


# set working directory
#setwd("D:\\Ramona_F\\PhD Projects\\Bats\\R Analysis\\DAPC") # uni
setwd("C:/Users/Jan_S/Desktop/2025-Mastomys")

# load library
library(adegenet)
library(ggplot2)
library(reshape2)
library(readr)


# load the matrix of z-values

x <- read_csv("SupertypeOutput/dist_matrix_20241113.csv")

x                           # must be symmetrical:


# Extract numeric columns and convert to a matrix
dist_matrix <- as.matrix(x[, -1])
# Set row names to the first column (if necessary)
rownames(dist_matrix) <- x[[1]]
# Symmetrize: average (i, j) and (j, i)
dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]
# Fill diagonal with 0
diag(dist_matrix) <- 0

x <- dist_matrix


### first lets check how stable the BIC-values are
## code from: https://grunwaldlab.github.io/Population_Genetics_in_R/clustering_plot.html

## set a number of clusters we want to check 

str(x)
head(x)

maxK <- 40
myMat <- matrix(nrow=40, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(x, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}


my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Groz-up", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)


# visualize

p2 <- ggplot(my_df, aes(x = K, y = BIC))
p2 <- p2 + geom_boxplot()
p2 <- p2 + theme_bw()
p2 <- p2 + xlab("Number of groups (K)")
p2




#####################################################################################################

############# AJ MHC II-Supertyping

#####################################################################################################

par(mfrow = c(2,2))

# when we create our final figures we need consistency, so we set a seed
# but when starting with DAPC it makes sense to also comment the set.seed() call out to see if results are very different with different seed
set.seed(1) 
grp <- find.clusters(x, max.n.clust = 50)

# Choose the number PCs to retain (>= 1): 
# 60
# Choose the number of clusters (>=2:  SUPERTYPES! 
# 19 ########### from MHC-Tools Result

set.seed(1)
dapc1 <- dapc(x, grp$grp)

# Choose the number PCs to retain (>=1):
# 15
# Choose the number discriminant functions to retain (>=1): 
# 10

# retaining too many PC-axis can lead to over-fitting and unstability in the membership probabilities, so start with 80 % of variance explained, then validate

# check if/how clustering can be improved
# a-score returns the optimal number of PC axis to keep

myascore<- optim.a.score(dapc1, n.sim=100) 

# optimal number of PCs: 9 (found in plot) 


# cross-validation
# another method for optimizing

xval <- xvalDapc(x, grp$grp, n.pca.max=40, training.set=0.9, 
                 result="groupMean", center=TRUE, scale=FALSE, 
                 n.pca=NULL, n.rep = 100, xval.plot = TRUE)
xval [2:6]
#$`Number of PCs Achieving Highest Mean Success` "10"
# Number of PCs Achieving Lowest MSE "10"

# number of PCs that should be retained: 35 # this seems too high when we compare with BIC-value curve


####### now repeat the dapc with the number of PCs suggested by the a-score

set.seed(1)
dapc1 <- dapc(x, grp$grp)

# Choose the number PCs to retain (>=1):
# 10
# Choose the number discriminant functions to retain (>=1): 
# 10

################## Also possible set own Supertype number (n.clust), i.e. result from MHCTools, but by putting in n.clust in line 90 also possible

#grp19 <- find.clusters( x , n.clust = 19)  

#dapc1 <- dapc(x, grp19$grp)


# plot results to see how clusters spread
par(mfrow = c(2,2))
scatter(dapc1, 1, 2)
scatter(dapc1, 1, 3)
scatter(dapc1, 2, 3)
scatter(dapc1, 1, 4)

summary(dapc1)


# check which alleles were assigned to which cluster / supertype
# this illustrates how well single alleles were assigned to a cluster
# ideally they should be mostly dark red (which means a high probability that this allele fits in this cluster)
par(mfrow = c(1,1))
assignplot(dapc1)
assignplot(dapc1, subset=1:50) # to look at a subset of the plot, e.g. the first 50 alleles
assignplot(dapc1, subset=51:90)
assignplot(dapc1, subset=91:131)
assignplot(dapc1, subset=1:100)
assignplot(dapc1, subset=101:200)
assignplot(dapc1, subset=200:297)

## this will plot allele membership in a different way and show membership probability to a supertype as a bar and coloured by supertype
compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:19),
          xlab="alleles", col=funky(19))



# but it is hard to tell which allele is which, so better use ggplot
# code annotated from: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

dapc.results <- as.data.frame(dapc1$posterior)
dapc.results$alleleID <- rownames(dapc.results)

dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Allele_ID","Assigned_Supertype","Posterior_membership_probability")

col=funky(19) # define colours 


p <- ggplot(dapc.results, aes(x= Allele_ID, y=Posterior_membership_probability, fill=Assigned_Supertype))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = col) 
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p <- p + theme(legend.position = "bottom")
p



# compoplot only for specific alleles ############################

#dapc.results <- as.data.frame(dapc1$posterior)
#dapc.results.panama <- dapc.results[1:133,] # select panama alleles
#dapc.results.panama$alleleID <- rownames(dapc.results.panama)


#dapc.results.panama <- melt(dapc.results.panama)

#colnames(dapc.results.panama) <- c("Allele_ID","Assigned_Supertype","Posterior_membership_probability")

#p1 <- ggplot(dapc.results.panama, aes(x= Allele_ID, y=Posterior_membership_probability, fill=Assigned_Supertype))
#p1 <- p1 + geom_bar(stat='identity') 
#p1 <- p1 + scale_fill_manual(values = col) 
#p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
#p1



#### now lets check which alleles are the most admixed ones #####################
## as admixed we define alleles with single-cluster membership probability < 85 % 

temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.85))) # we select here alleles that have no more than 85% of probability membership in a single cluster
temp


allele_names <-  names(temp)
names(temp)

# now lets plot only these admixed alleles

compoplot(dapc1, subset=temp, posi="bottomright",
          txt.leg=paste("Cluster", 1:19),
           col=funky(19), xlab = "alleles")
# Add custom x-axis labels with adjusted size and position
axis(1, 
     at = seq(1, length(temp) * 1.1, length.out = length(temp)), # Positions of the labels
     labels = allele_names, # Names of the alleles
     las = 2, # Rotate labels for better visibility
     cex.axis = 0.7, # Adjust label size (smaller than default)
     line = -.2) # Move labels downward (positive values shift them away)



###  Final assignment of alleles to clusters by posterior membership probability
dapc1$assign

# make into df
b<-dapc1$assign
c<-rownames(dapc1$posterior)

supertypes = data.frame(supertype=b, alleleID=c)
supertypes

# write csv of results for adding to the metadata
write.csv(supertypes,"19supertypes.long.csv") # finally write a csv file with alleles and assigned supertypes for our metadata!


