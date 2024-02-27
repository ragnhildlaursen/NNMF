# Function for batch updates of spatial nmf
#' @import Rcpp
#' @import RcppArmadillo

# sourceCpp("src/NMFspatial2.cpp")

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
#' @param location a matrix of the spatial locations of the observations 
#' @param lengthscale the lengthscale of the neighborhood. 
#' @param batch a vector of the same length as the observations to devide the observations into batches. Default is to have all observations in one batch
#' @param maxiter integer determining the maximum number of iterations
#' @param tolerance a small value to detmine the stopping criteria. default is 1e-8.
#' @param initial 
#' @param smallIter 
#' @param error_freq 
#' @param kernel_cutoff 
#' @param normalize 
#'
#' @return A list of the matrices derived by the factorization and the corresponding generalized Kullback-Leibler
#' \itemize{
#'  \item exposures - Non-negative matrix of dimension observations x noSignatures
#'  \item signatures - Non-negative matrix of dimension noSignatures x features, where rows sum to one.
#'  \item gkl - Smallest Value of the Generalized Kullback-Leibler
#'  \item gklvalues - The GKL values at each iteration from error_freq
#'  }
#' @export
#'
#' @examples
nnmf = function(data, noSignatures, location, lengthscale, batch = 1, maxiter = 10000, tolerance = 1e-8, initial = 5, smallIter = 100, error_freq = 10,kernel_cutoff = 0.5,normalize = TRUE){
    
    if(normalize){
        row_sum = rowSums(data)
        data = data/row_sum*mean(row_sum)
    }
    
    unique_batches = unique(batch)
    if(length(unique_batches) == 1){
        print("Everything is run in one batch")

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
    
    return(out)
}