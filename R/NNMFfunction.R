#' Compute Distance Between Locations
#' 
#' The distance between indexes of locations and all other locations in location matrix
#'
#' @param X matrix of location
#' @param index a vector of the indexes from which you want the distance to all other observations in the location matrix
#'
#' @return matrix of distances between the specified indexes and all other observations in the matrix of locations
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @import Matrix
dist_index = function(X,index){
  
  target <- X[index, ]  # Extract the index rows
  
  if(length(index) < 2){
    r <- apply(X, 1, function(xvar) sum((xvar - target)^2))
  }else{
    r <- apply(X, 1, function(x) apply(target, 1, function(y) sum((x - y)^2)))
    r = t(r)
  }
  
  if(any(r<0)){
    warning("Some distances were smaller than zero! Try scaling up the locations.")
    r[r<0] = 0
  }
  return(r)
}

#' The distance between two matrices of location
#'
#' @param X matrix of location
#' @param Y optional other matrix of location to take distance to
#'
#' @return symmetric matrix of distances.
#' @export
#'
dist_fun = function(X, Y = NULL){
  if(is.null(Y)){
    r = as.matrix(dist(X, diag = T, upper = T)^2)
  }else{
    if(ncol(X) != ncol(Y)){
      stop("The number of columns in X and Y need to be the same.")
    }
    
    X2 = rowSums(X^2)
    Y2 = rowSums(Y^2)
    X2 = matrix(X2, nrow = length(X2), ncol = length(Y2), byrow = F)
    Y2 = matrix(Y2, nrow = length(Y2), ncol = nrow(X2), byrow = F)
    
    r = X2 - 2*(X%*%t(Y)) + t(Y2)
  }
  
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
#' @param location a matrix of the spatial locations of the observations (observations x 2). If no location is specified regular NMF is used.
#' @param lengthscale the lengthscale of the neighborhood. If no lengthscale are specified then the lengthscale is set to 1/10 of the density of your points.
#' @param batch a vector of the same length as the observations to divide the observations into batches. Default is to have all observations in one batch.
#' @param maxiter integer determining the maximum number of iterations. Default is 10000
#' @param tolerance a small value to determine the stopping criteria. Default is 1e-8.
#' @param initial the number of intializations. Default is 1.
#' @param smallIter the number of iterations to run each initialization. Afterwards the initialization with the smallest error is run till convergence. Default is 100.
#' @param error_freq an interger determining how often to calculate the error. Default is 10.
#' @param kernel_cutoff a value between 0 and 1, where everything below this value is set to zero in the kernel. Default is 0.5.
#' @param normalize a logical value indicating whether the observations should be normalized. Default is TRUE.
#' @param not_sc a logical value that indicates if the data is not single cell data(default = TRUE). It is recommended to set to FALSE if the data is not single-cell.  
#'
#' @return A list of the matrices derived by the factorization and the corresponding generalized Kullback-Leibler
#' \itemize{
#'  \item weights - Non-negative matrix of dimension observations x noSignatures
#'  \item signatures - Non-negative matrix of dimension noSignatures x features, where rows sum to one.
#'  \item error - Final error of the Generalized Kullback-Leibler
#'  \item errorvalues - a vector of the length of maxiter. Includes the GKL values calculated at the frequency of error_freq
#'  \item avg_nn - The average number of nearest neighbors included in averaging.
#'  \item lengthscale - the lengthscale used for the analysis
#'  \item batch - The batches of the analysis, which is 1 if all were in one batch. 
#'  }
#' @export
#'
nnmf = function(data, noSignatures, location = NULL, lengthscale = NULL, batch = 1, maxiter = 1000, tolerance = 1e-10, initial = 3, smallIter = 50, error_freq = 10,kernel_cutoff = 0.1,normalize = TRUE, not_sc = FALSE){

  if(!is.matrix(data)){
    stop("The data needs to be of class matrix. \n")
  }
  if(sum(colSums(data) == 0) > 0){
    stop("Remove columns in the data that only contain zeroes. \n")
  }  
  
  if(sum(rowSums(data) == 0) > 0){
    stop("Remove rows in the data that only contain zeroes. \n")
  }  
    
  if(normalize){
    data <- sweep(data, 1, rowSums(data), FUN = "/") * ncol(data)
    cat('Normalized the data, so each row sums to ',ncol(data),'(the number of columns in the data). \n')
  }
    
    mean_nn = 0
    
    if(is.null(location)){
      cat("Running regular NMF, as no locations were specified. \n")
        out = nmfgen(data = data, noSignatures = noSignatures, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter, error_freq = error_freq)

    }else{
      
    if(nrow(data) != nrow(location)){
        stop("The number of rows in location must match the number of rows in data. \n")
      }
        
      unique_batches = unique(batch)
    
      if(length(unique_batches) == 1){
        cat("All ",nrow(data)," observations are run in one batch. \n")
        
        if(nrow(data) > 50000){
          stop("There is too many observation to run it in one batch. Use groupondist() to make batches with size 20000 or smaller. \n")
        }
        
        dist = dist_fun(location)
          
          if(is.null(lengthscale)){
            # est_ls = estimate_lengthscale(data = data, dist = dist, max_avg_nn = 25, column_ls = FALSE)
            # lengthscale = est_ls[which.min(est_ls[,2]),1]

            lengthscale = mean(apply(sqrt(dist),1,function(x) sort(x)[15]))
            
            cat("The length scale is set to", lengthscale, ". Adjust for a smaller or larger neighborhood after assessing results or use the function estimate_lengthscale() for a more appropriate length scale. \n")
          }
          
          # calculating covariance
          sigma = exp(-dist/(lengthscale^2))
          sigma[sigma < kernel_cutoff] = 0
          
          
          mean_nn = mean(rowSums(sigma > 0) - 1)
          sigma = sigma/rowSums(sigma)
          sigma <- Matrix(sigma, sparse = TRUE)
          weight = sigma
          
          out = nmfspatial(data = data, noSignatures = noSignatures, weight = weight, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter)

      }else{
        if(nrow(data) != length(batch)){
          stop("The length of batch must match the number of rows in data. \n")
        }
        
        weights = list()
        
        first_batch = TRUE

        batch_list <- list()
        
        for (i in unique_batches) {
          index <- which(batch == i)
          batch_list[[paste(i)]] <- index - 1  # Store numeric indices
          
          X = location[index,]
        
          dist = dist_fun(X)
            
          if(is.null(lengthscale)){
              lengthscale = mean(apply(sqrt(dist),1,function(x) sort(x)[15]))
              
              cat("The lengthscale is set to", lengthscale, ". Specify accordingly for a smaller or larger neighborhood after assessing results. \n")
          } 
            
          # calculating covariance
          sigma = exp(-dist/(lengthscale^2))
          sigma[sigma < kernel_cutoff] = 0
            
            
            if(first_batch){
              mean_nn = mean(rowSums(sigma > 0) - 1)
              first_batch = FALSE
            }else{
              mean_nn = 0.5*mean_nn + 0.5*mean(rowSums(sigma > 0) - 1)
            }
            sigma = sigma/rowSums(sigma)
            sigma <- Matrix(sigma, sparse = TRUE)
            weights[[paste(i)]] = sigma
        }

        if(initial == 1){
          print('error not happened')
          out = nmfspatialbatch2(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, error_freq = error_freq)
        }else{
          out = nmfspatialbatch(data = data, noSignatures = noSignatures, weight = weights, batch = batch_list, maxiter = maxiter, tolerance = tolerance, initial = initial, smallIter = smallIter, error_freq = error_freq)
          
      }

    }
    
    if(not_sc){
      out$exposures = out$exposures * ((data/(out$exposures%*%out$signatures)) %*% t(out$signatures))
    }
    }
    output = list()
    output$weights = out$exposures
    output$signatures = out$signatures
    output$error = out$gkl
    output$errorvalues = out$gklvalues
    output$avg_nn = mean_nn
    output$lengthscale = lengthscale
    output$batch = batch
    
    return(output)
}



