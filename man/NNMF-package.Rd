\name{NNMF}
\alias{NNMF-package}
\alias{NNMF}
\docType{package}
\title{
\packageTitle{NNMF}
}
\description{
\packageDescription{NNMF}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{NNMF}
\packageIndices{NNMF}
}
\author{
\packageAuthor{NNMF}

Maintainer: \packageMaintainer{NNMF}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the ~~
~~ R documentation directory ~~
\keyword{ package }
\seealso{
https://github.com/ragnhildlaursen/NNMF for more decription
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(NNMF)
# simulate data
data = matrix(rpois(1000*100,10), nrow = 1000) # load your count data here (observations x features)
location = matrix(runif(1000*2), nrow = 1000) # location of your observation (observations x 2)

res = nnmf(data = data, noSignatures = 2, location = location)

res$weights # The weights for each signature
res$signatures # The two signatures
}
