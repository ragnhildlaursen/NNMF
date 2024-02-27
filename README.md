# NNMF
Neighborhood Non-Negative Matrix Factorization

## Guide to use the package 
### Install
```
devtools::install_github("ragnhildlaursen/NNMF")
```

### Small user guide
```
library(NNMF)

data = matrix(rpois(10*20,10), nrow = 10) # load your count data here (observations x features)
location = matrix(runif(10*2), nrow = 10) # location of your observation (observations x 2)

res = nnmf(data = data, noSignatures = 2, location = location)

res$weights # The weights for each signature
res$signatures # The two signatures

```
