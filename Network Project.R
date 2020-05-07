
library(igraphdata)
library(igraph)
library(intergraph)

data("USairports")

#creating the .csv files
airport_verticies_df <- as_data_frame(USairports, what="vertices")
airport_edges_df <- as_data_frame(USairports, what="edges")

write.csv(airport_verticies_df, file="USairports Verticies.csv")
write.csv(airport_edges_df, file="USairports Edges.csv")

# I'm going to go with an undirected version
undir <- as.undirected(USairports, mode="each")

STundir <- asNetwork(undir) #creating the statnet version

detach(package:igraph)
library(statnet)
# Exploring

network.size(STundir) #755 airports
summary(STundir,print.adj=FALSE) #print.adj to FALSE prevents verbose output

#Vertex attributes
list.vertex.attributes(STundir)
# [1] "City"         "na"           "Position"     "vertex.names"
list.edge.attributes(STundir)
# [1] "Aircraft"   "Carrier"    "Departures" "Distance"   "na"         "Passengers" "Seats"

# network density
gden(STundir)
    # 0.08228082

# finding the diameter
lgc <- component.largest(STundir,result="graph")
gd <- geodist(lgc)
max(gd$gdist)
  # 8

# First plot
op <- par(mfrow=c(1,1), mar=c(1, 1, 1, 1))
gplot(STundir, vertex.col=2, vertex.cex=1, edge.col="light grey", edge.lwd=1,
      gmode="graph", main="US Airports")
par(op)


# COMMUNITY DETECTION
detach(package:statnet)
library(igraph)

# exploring modularity
table(E(undir)$Carrier)
E(undir)[1:10]$Departures
modularity(undir,(E(undir)$Carrier+1)) #error?
modularity(undir,(E(undir)$Departures+1)) #-0.008849193  not much
modularity(undir,(E(undir)$Seats+1)) # -0.008387755     not much
modularity(undir,(E(undir)$Passengers+1)) # -0.008935582  not much
modularity(undir,(E(undir)$Aircraft+1)) # 0.04757072  --> slight clustering here?
modularity(undir,(E(undir)$Distance+1)) # -0.004814224  not much

#trying many different community detecting algorithms and comparing.
cw <- cluster_walktrap(undir)
modularity(cw)
## 0.3548135
unique(membership(cw)) 
## very big output,many possible clusters (up to 112)
clp <- cluster_label_prop(undir)
modularity(clp)
## 0.2505687
unique(membership(clp)) #about 25 clusters
##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
cle <- cluster_leading_eigen(undir)
modularity(cle)
## 0.3884057
unique(membership(cle)) #about 10 clusters
## 1  7  8 10  2  9  3  4  5  6
cl <- cluster_louvain(undir)
modularity(cl)
## 0.4098616   --> best performing
unique(membership(cl)) #about 22 clusters 

# Let's plot what we found
op <- par(mfrow=c(1,1), mar=c(1, 1, 1, 1))
plot(cw, undir,vertex.color="blue",vertex.label=NA,vertex.size=3,edge.color="dark gray",
     main="Walktrap")
plot(clp, undir,vertex.color="blue",vertex.label=NA,vertex.size=3,edge.color="dark gray",
     main="Label Propagation")
plot(cle, undir,vertex.color="blue",vertex.label=NA,vertex.size=3,edge.color="dark gray",
     main="Leading Eigenvector")
plot(cl, undir,vertex.color="blue",vertex.label=NA,vertex.size=3,edge.color="dark gray",
     main="Louvain")
par(op)


# ERGM

detach(package:igraph)
library(statnet)
library(ergm)

# null model
DSmod0 <- ergm(STundir ~ edges, 
               control=control.ergm(seed=40)) # in this use of the 'control=' argument we set the seed.
summary(DSmod0)
# ==========================
#   Summary of model fit
# ==========================
#   
#   Formula:   STundir ~ edges
# 
# Iterations:  7 out of 20 
# 
# Monte Carlo MLE Results:
#       Estimate Std. Error MCMC % z value Pr(>|z|)    
# edges -4.10379    0.01483      0  -276.8   <1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Null Deviance: 394588  on 284635  degrees of freedom
# Residual Deviance:  47266  on 284634  degrees of freedom
# 
# AIC: 47268    BIC: 47278    (Smaller is better.)


DSmod1 <- ergm(STundir ~ edges + 
                 gwesp(0.7, fixed=TRUE),
               control=control.ergm(seed=40))
summary(DSmod1)
# ==========================
#   Summary of model fit
# ==========================
#   
#   Formula:   STundir ~ edges + gwesp(0.7, fixed = TRUE)
# 
# Iterations:  20 out of 20 
# 
# Monte Carlo MLE Results:
#                 Estimate Std. Error MCMC % z value Pr(>|z|)    
# edges           -13.9332     0.2974      1  -46.85   <1e-04 ***
# gwesp.fixed.0.7   5.3836     0.1400      1   38.46   <1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Null Deviance: 394588  on 284635  degrees of freedom
# Residual Deviance:  35709  on 284633  degrees of freedom
# 
# AIC: 35713    BIC: 35734    (Smaller is better.) 

#first look at goodness of fit
nullsim <- simulate(DSmod0, verbose = TRUE,seed = 5)
mod1sim <- simulate(DSmod1, verbose = TRUE, seed = 5)
rowgof <- rbind(summary(STundir ~ edges + degree(0:5) + triangle),
                summary(nullsim ~ edges + degree(0:5) + triangle),
                summary(mod1sim ~ edges + degree(0:5) + triangle) )
rownames(rowgof) <- c("USairports", "Null","Mod1")
rowgof
#             edges degree0 degree1 degree2 degree3 degree4 degree5 triangle
# USairports  4660      -5     126     111      83      49      44    26972
# Null        4524       0       0       0       2       9      26      671
# Mod1        3666     137      59      60      78      45      60    12178


