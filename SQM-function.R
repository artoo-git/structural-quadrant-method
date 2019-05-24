  ################################################################

# Structural Quadrant Method
#
# Botella, J. G. L. (2000). The structural quadrants method: 
# A new approach to the assessment of construct system complexity 
# via the repertory grid. Journal of Constructivist Psychology, 
# 13(1), 1-26.
# 
# ####################  CURRENTLY IS AT A DRAFT STAGE ##########
# I have written this function in the attempt to automatise the 
# calculation of the integration and differentiation indexes as 
# theorised by Galiffa and Botella in paper above. 
# the function returns the differentiation and integration indexes
# the method for analysis is obscure in the paper and I cannot
# replicate the results here. Therefore this function is currently
# in the status of DRAFT. 
#
# Diego Vitali

#################################################################

# The function here still has some defaults like rotate or nfactors 
# which I had put there for use with principal() (psych package). 
# They are unused with svd(), prcomp() or eigen()
#
# 

SQM <- function (x, min= "" , max = "", trim = 4, rotate = "none", nfactors = 2){

  X<-getRatingLayer(x, trim = trim)
  # calculate k as being the similarity parameter ( from Botella-Gallifa (2000) paper)
  k<- max-min
  # define number of rows and cols of the given n X m repertory grid
  rows<-dim(X)[1]
  cols<-dim(X)[2]


  ################################### SIMILARITY MATRICES (Botella-Gallifa 2000)
  
  # Given a N construct x M elements grid G: this extracts N similarity 
  # matrices S, one for each element, and constituted as N construct x 
  # (M-1) elements. Therefore, the matrix Sj for the element Ej is obtained 
  # by calculating the scoring similarity between Ej and every other element 
  # in every construct in G. 
  #
  #     sijk = r – |g ij – g ik | 
  #
  # r = max-min scoring 
  # gij = scoring of the element j - scoring of the element k in the construct i
  
  # the "result" variable is a matrix that will store the calculated similarity 
  # matrix in the for loop
  result<- matrix(0, nrow=rows, ncol=cols)
  # sim is an array ( actually a list) that will be filled with all the 
  # calculated similarity matrices
  sim<-list()
  for (e in seq(cols)){
    for (n in seq(rows)) {
      for (m in seq(cols)) {
        result[n,m] <- k - abs(X[n,e]-X[n,m]);
      }
    }
    sim[[paste('sim',e,sep='')]]<- result[,c(seq(cols)[seq(cols)!=e])]
    dimnames(sim[[e]]) <- list(rownames(X),colnames(X)[seq(colnames(X))!=e])
  }


  ################### INTEGRATION AND DIFFERENTIATION

  
  ### create integration and differentiation summary tables
  int<-matrix(0, nrow=rows, ncol=cols)
  dimnames(int) <- list(rownames(X),colnames(X))
  dif<-matrix(0, nrow=rows, ncol=cols)
  dimnames(dif) <- list(rownames(X),colnames(X))
  
  index<-length(sim)
  # Each similarity matrix Sj is now submitted to principal component analysis.
  for (s in seq(index)){
    # Selecting an element's similarity matrix in the list
    m<-as.matrix(sim[[s]])
    
    ############################################################################################
    ################ METHOD 1 USING EIGEN SPECTRAL DECOMPOSITION ###############################
    ############################################################################################
    #                                                                                          #
    # normalise similarity matrix. we need to do this if we are using the function principal().#
    # Principal() relies on eigen() which (I think) requires a normalised square symmetric     #
    # matrix:                                                                                  #
    # normalize() with option "1" normalise the matrix row-wise                                #
    #m <- normalize(m,normalize = 1)                                                           #
    #                                                                                          #
    # calculate correlation matrix of the transposed similarity matrix (we want the constructs)#
    # ( this is also necessary to avoid error due the similarity matrix being singular)         #
    #r<- cor(t(m))                                                                             #
    #                                                                                          #
    #                                                                                          #
    # Extracts factors for each similarity matrix Sj , and compute factor                      #
    # loadings for each construct c i:                                                         #
    #                                                                                          #
    # res<-principal(r, rotate = rotate, nfactors = nfactors)
    #                                                                                          #
    # ..extract loadings                                                                       #
    # loadings<- res$loadings[,1:2]   
    #                                                                                          #
    # ..extract eigenvalues                                                                    #
    # evalues<- res$values
    #                                                                                          #
    # *Or* we can use the standard eigen() and not principal() wich is a function included     #
    # in the psych package. The methods are equivalent.                                        #
    # res<-eigen(r)                                                                            #
    # loadings<- eigen(r)$vectors %*% sqrt(diag(eigen(r)$values, nrow = length(eigen(r)$values)))#
    #                                                                                          # 
    # ..and now eigenvalues. from eigen() are already eigenvalues and same is for principal()  #
    # evalues<- res$values                                                                     #
    #                                                                                          #
    ############################################################################################
    ############################################################################################
    
    
    ############################################################################################
    ############## METHOD 2 USING prcom() OR svd() ON RAW SIMILARITY MATRICES ##################
    ############################################################################################
    # EXPERIMENTAL: Gallifa says that they did factor analysed the raw similarity matrices     #
    # and not the respective correlation matrices. So I have tried use prcomp() with the raw   #
    # similarity matrices                                                                      #
    # normalize matrix row-wise 
    #                                                                                          #
    #m <- normalize(m,normalize = 1)                                                            #                                                          
    # transpose similarity so that it is factoranalysed using constructs as variables          #
    m<- t(m)
    #                                                                                          #
    res<-prcomp(m, center = T, scale. = F)                                                     #
    #                                                                                          #
    # get loadings                                                                             #
    loadings<-as.matrix(res$rotation)
    #                                                                                          #
    # calculate eigenvalues from sigma values                                                  #
    evalues<-res$sdev^2
    #
    ############################################################################################
    ############################################################################################
    
    
    #calculate percentage accunted by first and second factor 
    pvaff<-100*evalues[1]/sum(evalues)
    pvasf<-100*evalues[2]/sum(evalues)
    
    # Prepare the for loop and fill the differentiation and integration summaries
    l<-length(loadings[,1])
    for (n in seq(l)){
      a <-loadings[,1][n] * pvaff + loadings[,2][n] * pvasf
      b <-.3 * pvaff + .3 * pvasf
       ifelse(a<0 && abs(a)>b,
             dif[n,s]<-round(a,2),
             dif[n,s]<-as.numeric(NA))
      ifelse(a>b,
             int[n,s]<-round(a,2),
             int[n,s]<-as.numeric(NA))
    }
  }
  
  
  #########################  Need to write this part ################################
  # from Botella & Gallifa ( 2000): "Factor analyses do not indicate which pole of  #
  # ci is applied to ej [..] We consider that the same pole of ci is being applied  #
  # to ej and ek in G when |g ij – g ik | <= 1                                      #
  #                                                                                 #
  # (.. this I still need to write.... abs(X[n,e]-X[n,m])<=1) ...                   #
  ###################################################################################
  
  # N of differentiated elements (Del) and constructs (Dco)
  Del<-length(unique(which(dif < 0, arr.ind = TRUE)[,2]))
  Dco<-length(unique(which(dif < 0, arr.ind = TRUE)[,1]))
  
  # N of integrated elements (Iel) and N of differentiation constructs (Ico)
  Iel<-length(unique(which(int > 0, arr.ind = TRUE)[,2]))
  Ico<-length(unique(which(int > 0, arr.ind = TRUE)[,1]))
  
  # Calculate the integration and differentiation index
  indexInt<- Iel/cols * Ico/rows
  indexDiff<- Del/cols * Dco/rows
  
  res<- list("Integration index" = indexInt, "Integration summary" = int, 
             "Differentiation index" = indexDiff, "Differentiation summary" = dif)
  
  #########################################################################
  # Discussion                                                            #
  #                                                                       #
  # the summaries report all the respective diff and int weights          # 
  # but I am not sure they are right.. specially because the Authors      #
  # are extracting weights ne+2 whereas I am getting all weight below 100
  #
  # HOWEVER: using the second method (currently default) and so extracting
  # the factors applying principal component analysis directly on the 
  # transposed similarity matrices I seem to extract (almost) matching 
  # factors to those on the paper. 
  #                                                                       #
  #########################################################################
  
  return(res)
}

