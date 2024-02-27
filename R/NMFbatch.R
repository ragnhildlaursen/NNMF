# Function for batch updates of spatial nmf
#' @import Rcpp
#' @import RcppArmadillo

dist_fun = function(X){
    X2 = rowSums(X^2)
    X2 = matrix(X2, nrow = length(X2), ncol = length(X2), byrow = F)
    r = X2 - 2*X%*%t(X) + t(X2)
    if(any(r<0)){
        warning("Some distances were smaller than zero! Try scaling up the locations.")
        r[r<0] = 0
    }
        
    return(r)
}

dist_index = function(X,index){
    X2 = rowSums(X^2)
    r = sum(X[index,]^2) - 2*X[index,]%*%t(X) + X2
    if(any(r<0)){
        warning("Some distances were smaller than zero! Try scaling up the locations.")
        r[r<0] = 0
    } 
    return(r)
}


#' Creates batches of observations that are close in distance
#'
#' @param location a matrix of locations (observations x 2)
#' @param size an integer determining the size of the batches (overwrites no_group parameter)
#' @param no_groups an integer determining the number of batches 
#'
#' @return A vector of the same length as the number of locations
#' @export
#'
#' @examples
groupondist = function(location, size = NULL, no_groups = NULL){
    n = nrow(location)
    left = c(1:n)

    if(is.null(no_groups) & is.null(size)){
        stop("You must determine the size or number of groups")
    }

    if(is.null(size)){
        size = ceiling(n/no_groups)
    }

    i = 1
    batch_vec = rep(paste0("b",0),n)
    while(length(left) > size){ 
        start = sample(length(left),1)
        dist = dist_index(location[left,],start)
        batch_index = order(dist)[1:size]
        batch_vec[left[batch_index]] = paste0("b",i)
        i = i+1
        
        left = left[-batch_index]
    }
    return(batch_vec)
}


#' Neighborhood Non-Negative Matrix Factorization 
#' 
#' Finding factors that includes the information of the neighborhood through the spatial information.
#'
#' @param data a matrix of non-negative entries with the observations as rows and features as columns
#' @param noSignatures integer determining the rank of the factorization, which determines the number of signatures
#' @param location a matrix of the spatial locations of the observations (observations x 2) 
#' @param lengthscale the lengthscale of the neighborhood. 
#' @param batch a vector of the same length as the observations to divide the observations into batches. Default is to have all observations in one batch
#' @param maxiter integer determining the maximum number of iterations
#' @param tolerance a small value to determine the stopping criteria. default is 1e-8.
#' @param initial the number of intializations. default is 1.
#' @param smallIter the number of iterations to run each initialization. Afterwards the initialization with the smallest error is run till convergence.
#' @param error_freq an interger determining how often to calculate the error.
#' @param kernel_cutoff a value between 0 and 1, where everything below this value is set to zero in the kernel 
#' @param normalize a logical value indicating whether the observations should be normalized. 
#'
#' @return A list of the matrices derived by the factorization and the corresponding generalized Kullback-Leibler
#' \itemize{
#'  \item weights - Non-negative matrix of dimension observations x noSignatures
#'  \item signatures - Non-negative matrix of dimension noSignatures x features, where rows sum to one.
#'  \item error - Final error of the Generalized Kullback-Leibler
#'  \item errorvalues - a vector of the length of maxiter. Includes the GKL values calculated at the frequency of error_freq
#'  }
#' @export
#'
#' @examples
nnmf = function(data, noSignatures, location = NULL, lengthscale = NULL, batch = 1, maxiter = 10000, tolerance = 1e-8, initial = 1, smallIter = 100, error_freq = 10,kernel_cutoff = 0.5,normalize = TRUE){
    
    if(normalize){
        row_sum = rowSums(data)
        data = data/row_sum*mean(row_sum)
    }
    
  
    if(is.null(location)){
      cat("Running regular NMF, as no locations were specified.")
      out = nmfgen(data = data, noSignatures = noSignatures, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter, error_freq = error_freq)
    }else{
      
      if(nrow(data) != nrow(location)){
        stop("The number of rows in location must match the number of rows in data.")
      }
        
    if(is.null(lengthscale)){
      r1 = range(location[,1])
      r2 = range(location[,2])
      lengthscale = (r1[2] - r1[1])*(r2[2] - r2[1])/nrow(count)
      lengthscale = signif(lengthscale,1)/10
      cat("The lengthscale is set to", lengthscale, ". Specify accordingly for a smaller or larger neighborhood after assessing results.")
    }
    unique_batches = unique(batch)
    if(length(unique_batches) == 1){
      cat("All ",nrow(data)," observations are run in one batch.")
      if(nrow(data) > 50000){
        stop("There is too many observation to run it in one batch. Use groupondist() to make batches with size 20000 or smaller.")
      }

        dist = dist_fun(location)
        
        # calculating covariance
        sigma = exp(-dist/(lengthscale^2))
        sigma[sigma < kernel_cutoff] = 0

        weight = sigma/rowSums(sigma)

        out = nmfspatial(data = data, noSignatures = noSignatures, weight = weight, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter)

    }else{
        weights = list()
        batch_list = list()
        for(i in unique_batches){
            index = which(batch == i)
            batch_list[[i]] = index - 1

            X = location[index,]
        
            dist = dist_fun(X)
        
            # calculating covariance
            sigma = exp(-dist/(lengthscale^2))
            sigma[sigma < kernel_cutoff] = 0

            weights[[i]] = sigma/rowSums(sigma)
        }

        if(initial == 1){
            out = nmfspatialbatch2(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, error_freq = error_freq)
        }else{
            out = nmfspatialbatch(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter, error_freq = error_freq)
        }

    }
    }
    output = list()
    output$weights = out$exposures
    output$signatures = out$signatures
    output$error = out$gkl
    output$errorvalues = out$gklvalues
    return(output)
}



#' Find the top features for each signature
#'
#' @param signatures a matrix of the signatures (number of signatures x features)
#' @param feature_names a vector of names for the features.
#' @param ntop an integer of how many of the top features to output
#'
#' @return a matrix of top features for each signature. (number of signatures x ntop)
#' @export
#'
#' @examples
topfeatures = function(signatures, feature_names, ntop = 10){
  if(ncol(signatures) != length(feature_names)){
    stop("The feature_names need to be the same length as the number of columns in signatures.")
  }
  dat = t(signatures)
  dat_new=NULL
  for(ii in 1:nrow(dat)){
    rr=dat[ii,]
    m1=max(rr)
    m2=max(rr[-which(rr==m1)])
    mm=rep(m1, length(rr))
    mm[which(rr==m1)]=m2
    ns=rr*log((rr+1e-10)/(mm+1e-10))
    dat_new=rbind(dat_new, ns)
  }
  
  weight_topgene = NULL
  for(topic in 1:ncol(dat)){
    idx = order(dat_new[,topic], decreasing = T)
    weighting = feature_names[idx[1:ntop]]
    
    weight_topgene = rbind(weight_topgene,c(topic,weighting))
  }
  return(weight_topgene)
}
