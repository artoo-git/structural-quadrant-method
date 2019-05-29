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

# internal workhorse for SQM
#
# @param x               \code{repgrid} object.
# @param min             Minimum score should be define a priori (e.g for a range 1-7 
#                        min = 1). If not given the min is estimated on the basis of 
#                        the grid scoring (not recommended).
# @param max             Maximum score should be define a priori (e.g for a range 1-7 
#                        max = 7). If not given the max is estimated on the basis of 
#                        the grid scoring (not recommended).
# @param PCA.f           With \code{PCA.f="spectral"} SQM will call the function
#                        eigen() for spectral value decomposition. With \code{PCA.f="singular"}
#                        SQM will call the function prcomp() for singular value decomposition
# @param trim            The number of characters a construct (element) is trimmed to (default is
# @author                Diego Vitali
# @export
# @keywords internal
# @return                A list with four elements containing different steps of the 
#                        calculation.
#
# 

SQM <- function (x, min= "" , max = "", trim = 4, PCA.f= "spectral"){

  X<-getRatingLayer(x, trim = trim)
  
  sc <- getScale(x)
  if (is.na(min)){
    min <- floor(sc)
    cat("\nMinimum score not given. Estimated min based on grid scoring =", min)
  }
  if (is.na(max)){
    max <-  ceiling(sc)
    cat("\nMaximun score not given. Estimated Max based on grid scoring =", max)
  }
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
    

  if(PCA.f=="singular" | is.na(PCA.f)){ 
    
    ############################################################################################
    ############## METHOD 1 USING SINGULAR VALUE DECOMPOSITION VIA prcomp()   ##################
    ############################################################################################
    # EXPERIMENTAL: Gallifa says that they did factor analysed the raw similarity matrices     
    # and not the respective correlation matrices. So I have tried use prcomp() with the raw   
    # similarity matrices                                                                      
    # normalize matrix row-wise 

    # from "Principal Component Analysis in r:  an examination of the different functions and 
    # methods to perform PCA. Gregory B. Anderson:
    #
    # << The function prcomp() performs a principal component analysis using singular value decomposition 
    # of the data matrix.  Unlike the svd function, the data matrix does not need to be centered or 
    # scaled prior to analysis because these options can be specified as arguments in the prcomp 
    # command.  The arguments of the function include: formula (e.g., ~X1+X2), data(optional data 
    # rame containing the variables), subset (an optional vector used to select particular rows 
    # of the data matrix), na.action (a function to indicate how missing values should be treated), 
    # x (a matrix or dataframe to be used), retx(a logical statement to determine if the rotated 
    # variables should be returned), center (a logical statement to indicate that the data should
    # be centered at zero), scale (a logical statement to indicate  that the data should be 
    # scaled [default=FALSE]), tol(an argument to indicate that the magnitude below the provided 
    # value for components should be omitted), and newdata (an optional dataset to predict into).>>

    #m <- normalize(m,normalize = 1)                                                                                                                     
    # transpose similarity so that it is factoranalysed using constructs as variables          
    m<- t(m)
    #  
    # center: (the column means of the original data used to center each variable) 
    # scale:  (the original variance of each column used to scale each variable)                                                                                      #
    res<-prcomp(m, center = T, scale. = F)                                                     
    #                                                                                          
    # get loadings                                                                             
    loadings<-as.matrix(res$rotation)
    #                                                                                          
    # calculate eigenvalues from sigma values                                                  
    evalues<-res$sdev^2
    #
    ############################################################################################
    ############################################################################################
  
  } else if (PCA.f=="spectral"){

    ############################################################################################
    ################ METHOD 2 USING EIGEN SPECTRAL DECOMPOSITION ###############################
    ############################################################################################
    # 
    # from "Principal Component Analysis in r:  an examination of the different functions and 
    # methods to perform PCA. Gregory B. Anderson:
    #
    # << eigen() computes eigenvalues and eigenvectors using spectral decomposition of real or 
    # complex matrices.  The arguments for the function are as follows: 
    #    - X (a matrix to be used 
    #        [for PCA we specify either the correlation matrix or the covariance matrix])
    #    - symmetric (a logical argument; if TRUE the matrix is assumed to be symmetric 
    #        [if not specified, the matrix will be inspected for symmetry])
    #    - only.values (a logical argument; if TRUE only the eigenvalues are computed and returned). >>

    # We work on the similarity matrix:

    # we transpose similarity so that it is analysed using constructs as variables                                 #
    m<- t(m)                                                                                   
    #                                                                                          
    # calculate correlation matrix of the transposed similarity matrix (we want the constructs)
    # ( this is also necessary to avoid error due the similarity matrix being singular)         
    r<- cor(m)                                                                             
    #                                                                                          
    #                                                                                          
    res<-eigen(r)                                                                            
    loadings<- eigen(r)$vectors %*% sqrt(diag(eigen(r)$values, nrow = length(eigen(r)$values)))
    #                                                                                           
    # ..and now eigenvalues. from eigen() are already eigenvalues and same is for principal()  
    evalues<- res$values                                                                     
    #                                                                                          
    ############################################################################################
    ############################################################################################
  }  else{stop("\nthe PCA function selected is not recognised. quitting..")}

    
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

