#include <RcppArmadillo.h>
#include <cmath>        // std::abs
#include <tuple>
#include <iostream>


using namespace Rcpp;

double error(arma::colvec y, arma::colvec mu) {
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
      sum += pow(y[i] - mu[i],2);
  }
  return sum;
}

std::tuple<arma::mat, arma::mat, double> nmf1(arma::mat data, int noSignatures, int iter = 5000, double alpha = 10) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  arma::mat exposures(genomes, noSignatures,arma::fill::randu);
  
  arma::mat signatures(noSignatures, mutTypes, arma::fill::randu);
  
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));

    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );

    estimate = exposures * signatures;
    
    exposures = exposures % (((data/estimate) * arma::trans(signatures))/(((alphamat + data)/(alphamat + estimate))* arma::trans(signatures)));

    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;

    if(t - floor(t/5)*5 == 0){
      for(int r = 0; r<noSignatures; r++) {
        if(range(signatures.row(r)) <= 1e-10){
          Rcout << "collapsed. Try rerunning the function or check for zero columns or rows in your data.";
          break;
          // exposures.col(r) = arma::randg(genomes);
          // signatures.row(r) = arma::randg<arma::rowvec>(mutTypes);
        }
      }
    }
    
  }
  double gkl = error(arma::vectorise(data),arma::vectorise(estimate));
  
  return {exposures, signatures, gkl};
}

// [[Rcpp::export]]
List nmfgen_nb(arma::mat data, int noSignatures, int maxiter = 10000, double tolerance = 1e-8, int initial = 10, int smallIter = 100, int error_freq = 10, double alpha = 10) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  auto res = nmf1(data, noSignatures, smallIter);
  auto exposures = std::get<0>(res);
  auto signatures = std::get<1>(res);
  auto gklValue = std::get<2>(res);
  
  for(int i = 1; i < initial; i++){
    auto res = nmf1(data, noSignatures, smallIter);
    auto gklNew = std::get<2>(res);
    
    if(gklNew < gklValue){
      gklValue = gklNew;
      exposures = std::get<0>(res);
      signatures = std::get<1>(res);
      
    }
  }
  
  arma::mat estimate = exposures * signatures;
  
  double gklOld = error(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  arma::vec gklvalues(maxiter);
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    
    exposures = exposures % (((data/estimate) * arma::trans(signatures))/(((alphamat + data)/(alphamat + estimate))* arma::trans(signatures)));
    
    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    
    gklvalues.at(t) = gklOld;
    
    if(t - floor(t/error_freq)*error_freq == 0){
      gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
      
      if ((2*(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance) & (t > error_freq)){
        Rcout << "Total iterations:";
        Rcout << t;
        Rcout << "\n";
        break;
      }
      gklOld = gklNew;
    }
    
  }
  
  gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew,
                             Named("gklvalues") = gklvalues);
  return output;
}



// [[Rcpp::export]]
List nmfspatialbatch_nb(arma::mat data, int noSignatures, List weight, List batch, int maxiter = 10000, double tolerance = 1e-8, int initial = 10, int smallIter = 100, int error_freq = 10, double alpha = 10) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  int nobatches = batch.size();
  
  auto res = nmf1(data, noSignatures, smallIter);
  auto exposures = std::get<0>(res);
  auto signatures = std::get<1>(res);
  auto gklValue = std::get<2>(res);
  
  for(int i = 0; i < initial; i++){
    auto res = nmf1(data, noSignatures, smallIter);
    auto gklNew = std::get<2>(res);
    
    if(gklNew < gklValue){
      gklValue = gklNew;
      exposures = std::get<0>(res);
      signatures = std::get<1>(res);
      
    }
  }
  
  double gklOld = error(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  arma::vec gklvalues(maxiter);
  
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
  
    estimate = exposures * signatures;
    
    exposures = exposures % (((data/estimate) * arma::trans(signatures))/(((alphamat + data)/(alphamat + estimate))* arma::trans(signatures)));
   
    if(t - floor(t/2)*2 == 0){
    for(int b=0; b < nobatches; b++){
      arma::uvec batch_index = batch[b];
      arma::mat w_mat = weight[b];
      arma::mat exposures_batch = exposures.rows(batch_index);
      
      arma::colvec exp_sum = sum(exposures_batch,1);
      exposures_batch = exposures_batch.each_col() / exp_sum;
      exposures_batch = w_mat * exposures_batch;
      exposures_batch = exposures_batch.each_col() % exp_sum;

      exposures.rows(batch_index) = exposures_batch;
    }
    }
    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
 
    gklvalues.at(t) = gklOld;
    
    if(t - floor(t/error_freq)*error_freq == 0){
      gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
      if(t - floor(t/100)*100 == 0){
        Rcout << "Total iterations:";
        Rcout << t;
        Rcout << "  error:";
        Rcout << gklNew;
        Rcout << "\n";
      }
      if ((2*(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance) & (t > error_freq)){
        signatures = arma::normalise(signatures,1,1);
        Rcout << "Total iterations:";
        Rcout << t;
        Rcout << "\n";
        break;
      }
      gklOld = gklNew;
    }
    
  }
  
  gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew,
                             Named("gklvalues") = gklvalues);
  return output;

}

// [[Rcpp::export]]
List nmfspatialbatch2_nb(arma::mat data, int noSignatures, List weight, List batch, int maxiter = 10000, double tolerance = 1e-8, int error_freq = 10, double alpha = 10) {
  
  int nobatches = batch.size();
  
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;

  arma::mat exposures(genomes, noSignatures,arma::fill::randu);
  
  arma::mat signatures(noSignatures, mutTypes, arma::fill::randu);

  double gklOld = error(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  arma::vec gklvalues(maxiter);
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    
    exposures = exposures % (((data/estimate) * arma::trans(signatures))/(((alphamat + data)/(alphamat + estimate))* arma::trans(signatures)));
    
    if(t - floor(t/2)*2 == 0){
    for(int b=0; b < nobatches; b++){
      arma::uvec batch_index = batch[b];
      arma::mat w_mat = weight[b];
      arma::mat exposures_batch = exposures.rows(batch_index);

      arma::colvec exp_sum = sum(exposures_batch,1);
      exposures_batch = exposures_batch.each_col() / exp_sum;
      
      exposures_batch = w_mat * exposures_batch;
      
      exposures_batch = exposures_batch.each_col() % exp_sum;

      exposures.rows(batch_index) = exposures_batch;
    }
    }
    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
    
    gklvalues.at(t) = gklOld;
    
    if(t - floor(t/error_freq)*error_freq == 0){
      gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
      if(t - floor(t/100)*100 == 0){
      Rcout << "Total iterations:";
      Rcout << t;
      Rcout << "  error:";
      Rcout << gklNew;
      Rcout << "\n";
      }
      if ((2*(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance) & (t > error_freq)){
        Rcout << "Total iterations:";
        Rcout << t;
        Rcout << "\n";
        break;
      }
      gklOld = gklNew;
    }
    
  }
  
  gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew,
                             Named("gklvalues") = gklvalues);
  return output;

}

// [[Rcpp::export]]
List nmfspatial_nb(arma::mat data, int noSignatures, arma::mat weight, int maxiter = 10000, double tolerance = 1e-8, int initial = 5, int smallIter = 100, int error_freq = 10, double alpha = 10) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  auto res = nmf1(data, noSignatures, smallIter);
  auto exposures = std::get<0>(res);
  auto signatures = std::get<1>(res);
  auto gklValue = std::get<2>(res);
  
  for(int i = 1; i < initial; i++){
    auto res = nmf1(data, noSignatures, smallIter);
    auto gklNew = std::get<2>(res);
    
    if(gklNew < gklValue){
      gklValue = gklNew;
      exposures = std::get<0>(res);
      signatures = std::get<1>(res);
      
    }
  }
  
  double gklOld = error(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  arma::vec gklvalues(maxiter);
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
  
    estimate = exposures * signatures;
    
    exposures = exposures % (((data/estimate) * arma::trans(signatures))/(((alphamat + data)/(alphamat + estimate))* arma::trans(signatures)));
    
    if(t - floor(t/2)*2 == 0){
    arma::colvec exp_sum = sum(exposures,1);
    exposures = exposures.each_col() / exp_sum;
    exposures = weight * exposures;
    exposures = exposures.each_col() % exp_sum;
    }
    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
     
    gklvalues.at(t) = gklOld;
    if(t - floor(t/error_freq)*error_freq == 0){
    gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
    
    if ((2*(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance) & (t > error_freq)){
      Rcout << "Total iterations:";
      Rcout << t;
      Rcout << "\n";
      break;
    }
    gklOld = gklNew;
    }
  }
  
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew,
                             Named("gklvalues") = gklvalues);
  return output;
}

// [[Rcpp::export]]
List nmftrain_nb(arma::mat data, arma::mat exposures, arma::mat signatures, arma::mat sigma, arma::vec ls_vec, int maxiter = 5000, double tolerance = 1e-8, int error_freq = 10, double alpha = 10) {
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  int noSignatures = signatures.n_rows;

  arma::mat weight = sigma;
  
  double gklOld = error(arma::vectorise(data),arma::vectorise(estimate));
  double gklNew = 2*gklOld;
  arma::vec gklvalues(maxiter);
  
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  
  arma::mat estimate = exposures * signatures;
  
  for(int t = 0; t < iter; t++){
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    signatures = arma::normalise(signatures,1,1);
    
    signatures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
  
    estimate = exposures * signatures;
    
    exposures = exposures % (((data/estimate) * arma::trans(signatures))/(((alphamat + data)/(alphamat + estimate))* arma::trans(signatures)));
  
    if(t - floor(t/2)*2 == 0){
    arma::colvec exp_sum = sum(exposures,1);
    exposures = exposures.each_col() / exp_sum;
    for(int col = 0; col<noSignatures; col++) {
      if(ls_vec(col) > 0){
      weight = arma::pow(sigma,ls_vec(col));
      weight = arma::normalise(weight,1,1);
      exposures.col(col) = weight * exposures.col(col);
      }
    }
    exposures = exposures.each_col() % exp_sum;
    }
    exposures.transform( [](double val) {return (val < 1e-10) ? 1e-10 : val; } );
    
    estimate = exposures * signatures;
   
    gklvalues.at(t) = gklOld;
    if(t - floor(t/error_freq)*error_freq == 0){
      gklNew = error(arma::vectorise(data),arma::vectorise(estimate));
      
      if ((2*(gklOld - gklNew)/(0.1 + std::abs(2*gklNew)) < tolerance) & (t > error_freq)){
        Rcout << "Total iterations:";
        Rcout << t;
        Rcout << "\n";
        break;
      }
    gklOld = gklNew;
    }
  }
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = gklNew,
                             Named("gklvalues") = gklvalues);
  return output;
}