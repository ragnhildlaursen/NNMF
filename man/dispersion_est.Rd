% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dispersion_est}
\alias{dispersion_est}
\title{Estimation of the overdispersion parameter for Negative Binomial non-negative matrix factorization}
\usage{
dispersion_est(data, estimate, patient_specific = FALSE)
}
\arguments{
\item{data}{Numeric matrix of mutational counts data. Matrix size: no. of mutation types x no. of patients.}

\item{patient_specific}{Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.}

\item{k}{Number of signatures to be used for the non-negative matrix factorization}
}
\value{
Overdispersion parameter. Either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
}
\description{
Likehood estimation of the dispersion parameter in negative binomial using Newton-Raphson.
The overdispersion parameter can be either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
}
