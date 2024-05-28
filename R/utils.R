#' Find the top features for each signature
#'
#' @param signatures a matrix of the signatures (number of signatures x features)
#' @param feature_names a vector of names for the features.
#' @param ntop an integer of how many of the top features to output. Default is 10.
#'
#' @return a matrix of top features for each signature. (number of signatures x ntop)
#' @export
#'
#' @examples
topfeatures = function(signatures, feature_names, ntop = 10){
  if(ncol(signatures) != length(feature_names)){
    stop("The feature_names need to be the same length as the number of columns in signatures. \n")
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


#' Determine the average nearest neighbors of certain celltypes
#'
#' @param location A matrix of locations (observation x 2)
#' @param celltype A vector of the same length as the number of rows in the location matrix
#' @param nn        An integer determining the number of nearest neighbors
#' @param sampleid  A vector to split the location into different samples. Needs to be the same length as the number of rows in the location matrix or NULL (default) if all the observations come from the same sample. 
#'
#' @return
#' @export
#'
#' @examples
nn_adj = function(location,celltype,nn = 5,sampleid = NULL){
  sum_cell = table(celltype_total)
  
  dist_index = function(X,index){
    X2 = rowSums(X^2)
    if(length(index) < 2){
      r = sum(X[index,]^2) - 2*X[index,]%*%t(X) + X2
    }else{
      X2 = matrix(X2, nrow = length(index), ncol = length(X2), byrow = T)
      r = rowSums(X[index,]^2) - 2*X[index,]%*%t(X) + X2
    }
    
    if(any(r<0)){
      warning("Some distances were smaller than zero! Try scaling up the locations.")
      r[r<0] = 0
    } 
    return(r)
  }
  
  dist_cells = function(cell){
    if(is.null(sampleid)){
      r = dist_index(location, index = which(celltype == levels(celltype_total)[cell]))
      cells = apply(r,1,function(x) celltype[order(x)[2:(nn+1)]])
    }else{
      cells = c()
      for(i in 1:length(unique(sampleid))){
        location_sub = location[sampleid == unique(sampleid)[i],]
        celltype_sub = celltype[sampleid == unique(sampleid)[i]]
        
        idx = celltype_sub == levels(celltype_total)[cell]
        if(sum(idx) > 0){
          r = dist_index(location_sub, index = which(idx))
          cells_sub = apply(r,1,function(x) celltype_sub[order(x)[2:(nn+1)]])
          
          cells = cbind(cells,cells_sub)
        }
      }
    }
    
    nn = table(cells)
    nn_norm = 1/sum_cell
    nn_norm[names(nn)] = nn*nn_norm[names(nn)]
    nn_norm2 = nn_norm/sum(nn_norm)
    return(nn_norm2)
  }
  
  nn_adj_mat = sapply(1:length(levels(celltype)), dist_cells)
  colnames(nn_adj_mat) = paste0(levels(celltype),"-NN")
  return(nn_adj_mat)
}


#' Determine the length scale of your data - testing 10 values from 0 until a determined maximum value
#'
#' @param data        A matrix including your data 
#' @param location    A matrix with the same number of rows as the data
#' @param max_avg_nn  An integer determining the the average neighbor for the maximum length scale tested. The default is 20. 
#' @param max_pct     A number between 0 and 1 determining the percentage of data points to the maximum length scale 
#' @param dist        A symmetric distance matrix of the data points. Must have the same number of rows and columns as the number of rows in data. If this is specified the location is ignored. 
#' @param column_ls   A logical determining whether a length scale should be determined for each column in the data (TRUE) or default = FALSE, where only one length scale is determined.
#'
#' @return A data.frame of the length scales and the corresponding test error.
#' @export
#'
#' @examples
estimate_lengthscale = function(data, location = NULL, max_avg_nn = 20, max_pct = NULL, dist = NULL, column_ls = FALSE){

  
  if(is.null(dist)){
    if(is.null(location)){
      stop("You need to specify either location or dist (distance) matrix for your data points. ")
    }
    
    dist = dist_fun(location,location)
    diag(dist) = Inf
  }else{
    if(nrow(dist) != nrow(data) | ncol(dist) != nrow(data)){
      stop("The distance matrix need to have the same number of rows and columns as the data.")
    }
    
    diag(dist) = Inf
  }
  
  dist = sqrt(dist)
  min_val = min(dist)
  
  if(is.null(max_pct)){
      if(max_avg_nn<1 | round(max_avg_nn) != max_avg_nn){
        stop("max_avg_nn should be a positive integer.")
      }
      max_val = mean(apply(dist,1,function(x) sort(x)[max_avg_nn]))
  }else{
      if(max_pct > 1 | max_pct <= 0){
        stop("max_pct needs to be a value between 0 and 1.")
      }
      max_val = mean(apply(dist,1,function(x) sort(x)[floor(length(x)*max_pct)]))
  }
  
  lengthscale = seq(min_val, max_val, length.out = 10)
  
  data_norm = data/rowSums(data)
  
  test_error = c()
  for(j in 1:length(lengthscale)){
    if(lengthscale[j] > 0){
      sigma = exp(-dist^2/lengthscale[j]^2)
    }else{
      sigma = exp(-(dist + 1)*Inf)
    }
    
    sigma[sigma < 0.1] = 0
    
    r_sum = rowSums(sigma)
    sigma[r_sum == 0,] = 1
    
    weights = sigma/rowSums(sigma)
    if(column_ls){
      fit_group = colMeans((data_norm - weights%*%data_norm)^2)
      test_error = rbind(test_error,fit_group)
    }else{
      fit_group = mean((data_norm - weights%*%data_norm)^2)
      test_error = c(test_error,fit_group)
    }
    
  }
  
  return(data.frame(lengthscale, test_error, row.names=NULL))
}