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

SQM <- function (x, min= "" , max = "", trim = 4, PCA.m = "singular"){
  
  X<-getRatingLayer(x, trim = trim)
  
  sc <- getScale(x)
  
  if (is.null(min)){
    min <- floor(sc)
    cat("\nMinimum score not given. Estimated min based on grid scoring =", min)
  }
  if (is.null(max)){
    max <-  ceiling(sc)
    cat("\nMaximun score not given. Estimated Max based on grid scoring =", max)
  }
  # calculate k as being the similarity parameter ( from Botella-Gallifa (2000) paper)
  k<- max-min
  # define number of rows and cols of the given n X m repertory grid
  rows<-dim(X)[1]
  cols<-dim(X)[2]

  ####################################################################
  ########### SIMILARITY MATRICES (Botella-Gallifa 2000)
  ####################################################################
  # THIS SECTION SHOULD BE SOUND
  #
  # Given a N construct x M elements grid G: this extracts N similarity 
  # matrices S, one for each element, and constituted as N construct x 
  # (M-1) elements. Therefore, the matrix Sj for the element Ej is obtained 
  # by calculating the scoring similarity between Ej and every other element 
  # in every construct in G. 
  #
  #     s_ijk = r - |g_ij - g_ik | 
  #
  # r = max-min scoring 
  # g_ij = scoring of the element j - scoring of the element k in the construct i
  
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


  ########################################################################
  ######### INTEGRATION AND DIFFERENTIATION INDEXES
  ########################################################################
  
  ### create empty integration and differentiation summary tables

  int<-matrix(0, nrow=rows, ncol=cols)
  dimnames(int) <- list(rownames(X),colnames(X))
  dif<-matrix(0, nrow=rows, ncol=cols)
  dimnames(dif) <- list(rownames(X),colnames(X))
  
  index<-length(sim)

  # Each similarity matrix Sj in the array is now submitted to principal 
  # component analysis.
  
  for (s in seq(index)){
  
    # Selecting an element's similarity matrix in the list
  
    m<-as.matrix(sim[[s]])


    # BOTH THE CHOICES BELOW AND THE CONFIGURATION OF EACH PCA METHOD ARE "SKETCHY" 
    # .. In fact some aspects I don't think make much sense they way they are described in the paper. 
    
    # (cf. p. 10-11, Gallifa & Botella 2000)
    # << Each similarity matrix is now submitted to factor analysis. Among all possible algorithms 
    #    that extract the matrix factors, we used the one by H orst (1965, section 12.5, pp. 624–625), 
    #    which is also the one used by Bell (1987) in his computer program for the analysis of repertory 
    #    grids G-Pack (version 3.0) >>
    #
    #    please note that both methods mentioned by the authors are are eigen value decomposition (PCA).
    #
    # << We extract two factors for each similarity matrix s_j , and compute factor loadings 
    #    for each construct c_i . Factor loadings for every c_i indicate the level of similarity 
    #    between the scorings of e_j and the rest of elements in G (standard rep. grid) considered 
    #    as a whole. Positive factor loadings indicate that the scoring of e_j in c_i is similar 
    #    to the scoring of all the other elements in c_i . In this case, the use of c_i when rating 
    #    e_j is very similar to its use when rating all the other elements in G. (Note, however, 
    #    that factor analyses do not indicate which pole of c_i is applied to e_j ; this information 
    #    needs to be obtained from the original ratings in G, as we will discuss later in more detail). 
    #    Negative factor loadings indicate that the scoring of e_j in c_i is dissimilar to the scoring 
    #    of all the other elements in c_i>>
    #

  if(PCA.m =="singular" | is.na(PCA.m)){ 

    ############################################################################################
    ############## METHOD 1 USING SINGULAR VALUE DECOMPOSITION VIA prcomp() 
    ############################################################################################
    # EXPERIMENTAL: Gallifa says that they did factor analysed the raw similarity matrices     
    # and not the respective correlation matrices. So I have tried use prcomp() with the raw   
    # similarity matrices                                                                      
    
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

    #m <- normalize(m,normalize = 1)  # I don't think this is needed                                                                                                               
    
    # transpose similarity so that it is submitted to PCA using constructs as variables
    m<- t(m)
    m<-cor(m)
    # center: (the column means of the original data used to center each variable) 
    # scale:  (the original variance of each column used to scale each variable)
    
    res<-prcomp(m, center = T, scale. = F)                                                     
    #
    # calculate eigenvalues from sigma values                                                  
    evalues<-res$sdev * res$sdev
    # pull eigenvectors 
    evectors<-as.matrix(res$rotation)
    #
    # Now, Loadings should be attained by: Eigenvectors * sqrt(Eigenvalues), but it maybe the case that
    # here Gallifa and botella only meant just the coeficient of the transformation (vectors)
    #
    # -------> Which one doe we use? <----------
    #
    # actual loadings:
    #loadings<- evectors %*% sqrt(diag(evalues, nrow = length(evalues)))
    #
    # or eigenvectors?                                                                          
    loadings<-evectors
    #       
    #
    ############################################################################################
    ############################################################################################
  
  } else if (PCA.m=="spectral"){

    ############################################################################################
    ################ METHOD 2 USING EIGEN SPECTRAL DECOMPOSITION
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
    r<- cor(m)
    #                                                                                          
    #
    # calculate eigenvalues.   
    evalues<- eigen(r)$values
    
    # calculate eigenvectors. 
    evectors<- eigen(r)$vectors                                                                   
    #                           
    #                           
    # Now, Loadings should be attained by: Eigenvectors * sqrt(Eigenvalues), but it maybe the case that
    # here Gallifa and botella only meant just the coeficient of the transformation (vectors)
    #
    # -------> Which one doe we use? <----------
    #
    # actual loadings:
    # loadings<- evectors %*% sqrt(diag(evalues, nrow = length(evalues)))
    #
    # or eigenvectors ?
    loadings<-evectors
    #                                                                              
    ############################################################################################
    ############################################################################################
    
  }  else if (PCA.m == "EFA"){
    #
    ############################################################################################
    ################ METHOD 3 USING FACTOR ANALYSIS  (discouraged)
    ############################################################################################
    #
    # I do not think this is the kind of analysis that gallifa and botella meant this method should
    # be removed
    library("psych")
    # we transpose similarity so that it is analysed using constructs as variables                                 #
    m<- t(m)
    #
    res<-fa(r = cor(m), nfactors = 2)
    loadings<-res$loadings
    #
    ############################################################################################
    ############################################################################################
    
    
  }else{stop(" The PCA function selected is not recognised. Please use 'singular' for singular value
          decomposition or 'spectral' for spectral value decomposition. quitting..")}

  ############################################################################################
  ################ WRAPPING UP RESULTS AND BUILDING TABLES 
  ############################################################################################
    
  # THIS SHOULD BE SOUND AS WELL (provided the pca is fine as it is)

    #calculate percentage accunted by first and second factor 
    # This is only different if Exploratory factor analysis is used 
    if (PCA.m == "EFA"){
      pvaff<-res$Vaccounted[4,1]
      pvasf<-res$Vaccounted[4,2]
       
    }else{
      pvaff<-100*evalues[1]/sum(evalues)
      pvasf<-100*evalues[2]/sum(evalues)
    }
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
  
  
  #########################  This part is rather Sketchy and curious ################################
  # from Botella & Gallifa ( 2000): "Factor analyses do not indicate which pole of  #
  # ci is applied to ej [..] We consider that the same pole of ci is being applied  #
  # to ej and ek in G when |g ij – g ik | <= 1                                      #
  #                                                                                 #
  # (.. this I still need to write.... abs(X[n,e]-X[n,m])<=1) ...  
  #
  # However .. I cannot really make sense of this passage... the signs of the loading
  # on components extracted only indicate their orthogonality relatively to the matrices
  # analysed and they don't have any other mathematical "meaning".
  #
  ###################################################################################
  
  ############################################################################################
  ################ PRINTING OFF 
  ############################################################################################
    # 
  # N of differentiated elements (Del) and constructs (Dco)
  Del<-length(unique(which(dif < 0, arr.ind = TRUE)[,2]))
  Dco<-length(unique(which(dif < 0, arr.ind = TRUE)[,1]))
  
  # N of integrated elements (Iel) and N of differentiation constructs (Ico)
  Iel<-length(unique(which(int > 0, arr.ind = TRUE)[,2]))
  Ico<-length(unique(which(int > 0, arr.ind = TRUE)[,1]))
  
  # Calculate the integration and differentiation index
  indexInt<- Iel/cols * Ico/rows
  indexDiff<- Del/cols * Dco/rows
  
  res<- list("Differentiation index" = indexDiff, "Differentiation summary" = dif,
             "Integration index" = indexInt, "Integration summary" = int)
  
  #########################################################################
  # Discussion                                                            
  #                                                                       
  # the summaries report all the respective diff and int weights           
  # but I am not sure they are right.. specially because the Authors      
  # are extracting weights ne+2 whereas I am getting all weight below 100
  #
  # HOWEVER: using spectral value decomposition and so extracting
  # the components with PCA directly on the 
  # transposed similarity matrices seem to yield  results similar to the
  # but I need to review the methods for PCA and ensure they make sense as
  # they are                                                               
  #########################################################################
  
  return(res)
}

