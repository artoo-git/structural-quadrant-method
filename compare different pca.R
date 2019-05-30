
sim<-function(n.var=5,n.obs=1000,n.reps=1000){
  results<-data.frame(matrix(numeric(0),n.reps,3))
  colnames(results)<-c("Values","Vectors","Scores")
  for(h in 1:n.reps){
    #Function to create data
    rPCA<-function(n.var,n.obs){
      mu<-rnorm(n.var)
      R<-rcorrmatrix(n.var,alphad=1)
      X<-rmvnorm(n=n.obs,mean=mu,sigma=R)
      return(X)
    }
    data<-rPCA(n.var=n.var,n.obs=n.obs)
    ##################################
    #Create Objects from PCA analysis
    ##################################
    #Analysis using metrics calculated from the function eigen
    pc1<-PC(data,scale=T,rm.na=T,print.results=F)
    #Analysis using the R function prcomp 
    pc2<-prcomp(data,center=TRUE,scale=TRUE) 
    #Analysis using the R function princomp
    pc3<-princomp(data,cor=TRUE) 
    #Analysis using the function PCA from the package FactoMinR
    pc4<-PCA(data, scale.unit=TRUE, ncp=ncol(data), graph=F)
    #Analyis using the function pca from the package pcaMethods
    pc5<-pca(data,method="svd",scale="uv",nPcs=ncol(data))  
    ############################################
    #Compare Eigenvalues from different methods
    ############################################
    eigenvalues<-data.frame(pc1$values,pc2$sdev^2,pc3$sdev^2,pc4$eig$eigenvalue,pc5@sDev^2)
    colnames(eigenvalues)<-c("PC","prcomp","princomp","PCA","pca")
    rownames(eigenvalues)<-colnames(data)
    valuediff<-matrix(numeric(0),nrow(eigenvalues),ncol(eigenvalues))
    for(i in 1:ncol(data)){
      for(j in 1:ncol(eigenvalues)){
        valuediff[i,j]<-eigenvalues[i,j]-rowMeans(eigenvalues)[i]
      }
    }
    results$Values[h]<-mean(valuediff) #Should be close to Zero
    #############################################
    #Compare Eigenvectors from different methods#
    #############################################
    eigenvectors<-list(PC=pc1$vectors,prcomp=pc2$rotation,princomp=pc3$loadings,PCA=pc4$svd$V,pca=pc5@loadings)
    ev<-array(0, dim=c(5,ncol(data),ncol(data)))
    ev[1,,]<-eigenvectors$PC
    ev[2,,]<-eigenvectors$prcompev[3,,]<-eigenvectors$princomp
    ev[4,,]<-eigenvectors$PCAev[5,,]<-eigenvectors$pca 
    vectordiff<-array(0, dim=c(5,ncol(data),ncol(data)))
    for(i in 1:5){
      for(j in 1:ncol(data)){
        for(k in 1:ncol(data)){vectordiff[i,j,k]<-abs(ev[i,j,k])-mean(abs(ev[,j,k]))
        }
      }
    }
    results$Vectors[h]<-mean(vectordiff) #Should be close to Zero
    #######################################
    #Compare Scores from different methods#
    #######################################
    scores<-list(PC=pc1$scores,prcomp=pc2$x,princomp=pc3$scores,PCA=pc4$ind$coord,pca=pc5@scores)
    sc<-array(0,dim=c(5,nrow(data),ncol(data)))
    sc[1,,]<-scores$PC
    sc[2,,]<-scores$prcomp
    sc[3,,]<-scores$princomp
    sc[4,,]<-scores$PCA
    sc[5,,]<-scores$pcascorediff<-array(0,dim=c(5,nrow(data),ncol(data)))
    for(i in 1:5){
      for(j in 1:nrow(data)){
        for(k in 1:ncol(data)){scorediff[i,j,k]<-abs(sc[i,j,k])-mean(abs(sc[,j,k]))
        }
      }
    }
    results$Scores[h]<-mean(scorediff)
  }
  return(results)
}

