# NNMF
Neighborhood Non-Negative Matrix Factorization

## Guide to use the package 
So far the package contain three different functions. 
 - *nnmf* which include the neighborhood nmf method to run your data together with the location.
 - *groupondist* creates batches of your data from the locations. It makes batches of observations that is close in distance. Recommended to use when the dataset has more than 20.000 observations.
 - *topfeatures* returns the top important features for your signatures.
### Install
```
devtools::install_github("ragnhildlaursen/NNMF")
```

### Small user guide for simulated data
```
library(NNMF)
# simulate data
data = matrix(rpois(1000*100,10), nrow = 1000) # load your count data here (observations x features)
location = matrix(runif(1000*2), nrow = 1000) # location of your observation (observations x 2)
```
If you have a very large dataset (more than 20.000 observations) it is recommended to split your dataset into batches. Use the groupondist function to do this.
```
batch_id = groupondist(location,size = 100)
```
Now you can run your model with the specified batches and you can start by using the predined lengthscale, which you can later try to change after assessing your results.
```
res = nnmf(data = data, noSignatures = 2, location = location, batch = batch_id)

res$weights # The weights for each signature
res$signatures # The two signatures

```
