######################Required packages#####################

require(Matrix)
require(pheatmap)

######################Parametros##############################
#Quantity of clusters, variables and instances
nc<-3
nv<-500
np<-300
# Distribution of instances in clusters and noise parameters
dp<-list(1:150,160:200,210:300)
noise_cluster<-c(0.5,0,0)
noise_general<-0.1
length(unlist(dp))
length(dp)

#pmatrix<-matrix(c(rep(0,nc-1),rnorm((nc-1)*(nc-1))),nc,(nc-1),byrow = TRUE)
#pmatrix  
#rankMatrix(pmatrix)
#clusterdistances<-as.matrix(dist(pmatrix,upper=TRUE,diag = TRUE))


# Distance matrix between clusters

clusterdistances<-matrix(c(0,1,1,1,0,2,1,2,0),3,3)
clusterdistances

######################Functions###############################
fromdtos<-function(B){
  M<-matrix(0,nrow(B),ncol(B))
  for(i in 1:nrow(B)){
    for(j in 1:ncol(B)){
      M[i,j]<-(B[i,j]-(B[1,i]+B[1,j]) )*(-0.5)
    }
  }
  return(as.matrix(M))
}  
normalizar<-function(v){
  return(v/(sqrt(sum(v*v))))
}
fromdistancetoset<-function(DIS){
  nd<-ncol(DIS)-1
  B<-DIS*DIS
  diagonalizacion<-eigen(fromdtos(B))
  if(sum(diagonalizacion$values<0)==0){
    vectores<-as.matrix(apply(diagonalizacion$vectors,2,normalizar))
    valores<-as.numeric(diagonalizacion$values)
    valores<-as.numeric(sort(valores,decreasing = TRUE))
    valores_sel<-as.numeric(valores[1:nd])
    P<-as.matrix(vectores[,order(valores,decreasing = TRUE)[1:nd]],nrow=nrow(B),ncol=nd)
      if(nd>1){
        D<-as.matrix(sqrt(diag(valores_sel)),nd,nd)
      }else{
        D<-valores_sel[1]
      }  
    proy<-P%*%D
    proy[1:(nd+1),1]<-proy[1:(nd+1),1]+1
    return(as.matrix(proy))
      
  }else{
    print("Error: invalid distance matrix")
    break
  }
}
builtbasis<-function(nc,nv){
  res<-nv%/%nc
  x<-seq(from=-pi,to=pi,length.out = res)
  for(i in 1:nc){
    if(i==1){
      z<-c(z,(1+cos(x))/2,rep(0,(nv-res)))
    }else if(i==nc){
      z<-c(z,rep(0,(nv-res)),(1+cos(x))/2)
    }else{
      z<-c(z,rep(0,res*(i-1)),(1+cos(x))/2,rep(0,((res*(nc-i))+(nv%%nc))))
    }
  }
  z<-z*(sqrt(2/(3*pi)))
  basis<-matrix(z,ncol=nv,nrow = nc,byrow = TRUE)
  return(as.matrix(basis))
}

#####################Building############################

coordinates<-fromdistancetoset(clusterdistances)
basis<-builtbasis(nc-1,nv)

data<-matrix(0,np,nv)

for(i in 1:np){
  m<-c()
  for(j in 1:length(dp)){
    if( i %in% dp[[j]]){
      m<-c(m,j)
    }
  }
  if(length(m)>1){
    coordinates_noisy<-apply(coordinates[m,],2,mean)+rnorm(mean=0,sd=0.1*mean(noise_cluster[m]),n=(nc-1))
  }else if(length(m)==1){
    coordinates_noisy<-coordinates[m,]+rnorm(mean=0,sd=0.1*noise_cluster[m],n=(nc-1))
  }
  if(length(m)==0){
    coordinates_noisy<-rnorm(mean=0,sd=5*noise_general,n=(nc-1))
  }
  data[i,]<-coordinates_noisy%*%basis
}
data<-data+matrix(rnorm(mean=0,sd=noise_general,n=np*nv),np,nv)

####################Display and export################
data_distance<-as.matrix(dist(data,upper=TRUE,diag=TRUE))
pheatmap(data,color=colorRampPalette(c("yellow","red","black"))(1000),cluster_cols = FALSE, cluster_rows = FALSE)
pheatmap(data_distance,color=colorRampPalette(c("yellow","red","black"))(1000),cluster_cols = FALSE, cluster_rows = FALSE)
