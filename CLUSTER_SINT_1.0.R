

require("pheatmap")
rm(list=ls())
nombre<-"X.csv"

# Parameters introduction stage


# General Parameters
n<-1500
m<-800
noisefactorpat_1<-0.15
noisefactorpat_2<-0.15
noisefactor<-1
margingenes<-0
marginpacientes<-0
# Cluster 1 parameters
mu<-0
pc1<-1:400
ng1<-length(pc1)
p<-400
alpha<-0
#Cluster 2 parameters
delta<-0
pc2<-41:80
ng2<-length(pc2)
q<-40
beta<-0


############################################################################

# Computing and Building datasets stage


# Reescaling of noise parameters.
noisefactorpat_1<-noisefactorpat_1*6
noisefactorpat_2<-noisefactorpat_2*6
noisefactor<-noisefactor*6


# A matrix of dimensions is created first
# required populated by samples of a normal distribution
# of mean 0 and variance equal to noisefactor.

noisematrix<-matrix(rnorm((m+2*marginpacientes)*(n+2*margingenes),mean=0,sd=noisefactor),m+2*marginpacientes,n+2*margingenes)


#h_k:[0,2*pi]---->[0,1] k e {0,...,n1-1}


# The function h_k acts as a real function of the type
# ((1 + cos (xK + M) / 2) ^ (1 / p) modulated by the parameters of the
# cluster 1. h_k assigns values
# high to the first variables of the individuals
# belonging to the first cluster and low or null values
# the last variables.

h_k<-function(s){
  alpha_rand<-rnorm(1,mean=alpha,sd=0.01)
  if((s>=0)&&(s<((1+mu)*pi*(1+alpha_rand)))){
    return(((1+cos(((s/(1+mu))-alpha_rand*pi)))/2)^{1/p})
  } else if((s>=((1+mu)*pi*(1+alpha_rand)))&&(s<=2*pi)){
    return(0)
  }else{
    return("error dominio h_k")
  }
}
#g_k:[0,2*pi]---->[0,1] k e {0,...,n2-1}


# g_k is the fully analogous version to h_k that provides the
# patterns of the second cluster. The modulation is the same as
# for h_k. g_k assigns high values to the last variables
# of the individuals of the second cluster and low values or
# null to the first variables.

g_k<-function(s){
  beta_rand<-rnorm(1,mean=beta,sd=0.01)
  if((s<(1-delta)*pi*(1-beta_rand))&&(s>=0)){
    return(0)
    
  } else if(s>=(1-delta)*pi*(1-beta_rand)&&(s<=2*pi)){
    return(((1+cos(((2*pi-s)/(1+delta))-beta_rand*pi))/2)^(1/q))
    
  }else{
      return("error dominio g_k")
  }
}


# matrixcdef is memory space dedicated to the construction of the
# final matrix of synthetic data. It is initialized now as an array
# null of the size agreed by patients and variables that contain
# desired patterns and with marginal and marginal parameters,
# which are the number respectively of individuals and variables that do not
# present no pattern, in order to present more realism.
# Next we add noisematrix to then include the
# patterns where appropriate, making use of h_k and g_k.

matrixcdef<-matrix(0,m+2*marginpacientes,n+2*margingenes)+noisematrix

x<-seq(from=0,to=2*pi,length.out = n)
if(ng1>1){
cluster1<-function(){
  c1m<-matrix(0,ng1,n)
  for(i in 0:(ng1-1)){
    v<-sapply(x, h_k )
    c1m[i+1,]<-v
  }
  rf<-sample(ng1)
  c1m<-c1m[rf,]
  return(c1m)
}
p1<-cluster1()

# Patterns are altered using modulated noise matrices in this
# case by noisefactorpat_1 and by noisefactorpat_2 in the case of the second
#cluster.

noisematrixcluster1<-matrix(rnorm(ng1*n,mean=0,sd=noisefactorpat_1),ng1,n)
p1<-p1+noisematrixcluster1
#pheatmap(p1,color=colorRampPalette(c("grey","blue"))(100),cluster_cols = FALSE, cluster_rows = FALSE)
matrixcdef[marginpacientes+pc1,((margingenes+1):(n+margingenes))]<-p1
}

if(ng2>1){
cluster2<-function(){
  c2m<-matrix(0,ng2,n)
  for(i in 0:(ng2-1)){
    v<-sapply(x, g_k)
    c2m[i+1,]<-v
  }
  rf<-sample(ng2)
  c2m<-c2m[rf,]
  return(c2m)
}
p2<-cluster2()
noisematrixcluster2<-matrix(rnorm(ng2*n,mean=0,sd=noisefactorpat_2),ng2,n)
p2<-p2+noisematrixcluster2
matrixcdef[(marginpacientes+pc2),((margingenes+1):(margingenes+n))]<-p2
}


# In the event that the list of individuals desired in the
# cluster 1 and in 2 have non-zero intersection the individuals of said
# intersection will be endowed with a pattern corresponding to the average of
# the characteristic patterns of each cluster, simulating a
# real intersection between the two subgroups.

pci<-intersect(pc1,pc2)
ngi<-length(pci)
if(ngi!=0){
  cim<-matrix(0,ngi,n)
  for(i in 0:(ngi-1)){
    random<-runif(1,0.25,0.75)
    v<-(random)*sapply(x, h_k)+(1-random)*sapply(x, g_k)
    cim[i+1,]<-v
  }
  noisematrixclusteri<-matrix(rnorm(ngi*n, mean=0, sd=(noisefactorpat_1+noisefactorpat_2)/2),ngi,n)
  pi<-cim+noisematrixclusteri
  matrixcdef[(marginpacientes+pci),((margingenes+1):(margingenes+n))]<-pi
}
# matrixf is finally returned, the forced version to be data frame of
# matrixcdef (synthetic data matrix)

matrixf<-as.data.frame(matrixcdef)
#Heatmap of the data
pheatmap(matrixf,color=colorRampPalette(c("black","red","yellow"))(1000),cluster_cols = FALSE, cluster_rows = FALSE)
#Heatmap of the correlation
pheatmap(as.matrix(dist(matrixf)),color=colorRampPalette(c("yellow","red","black"))(1000),cluster_cols = FALSE, cluster_rows = FALSE)

print("Guardar?")
p<-0
if(p==1){
  write.csv(matrixf,nombre,row.names=FALSE)
}


require(PINSPlus)

datalist<-list(as.matrix(matrixf),as.matrix(matrixf))

result<-SubtypingOmicsData(datalist)
 
