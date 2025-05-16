# NNMF
Neighborhood Non-Negative Matrix Factorization

## Guide to use the package 
So far the package contain three different functions. 
 - *nnmf* which include the neighborhood nmf method to run your data together with the location.
 - *groupondist* creates batches of your data from the locations. It makes batches of observations that is close in distance. Recommended to use when the dataset has more than 20.000 observations.
 - *topfeatures* returns the top important features for your signatures.
## Install
```
devtools::install_github("ragnhildlaursen/NNMF")
```

## Small user guide for simulated data
Below we load the package and simulate a dataset of size 1000x100. 
```
library(NNMF)
# simulate data
data = matrix(rpois(1000*100,10), nrow = 1000) # load your count data here (observations x features)
location = matrix(runif(1000*2), nrow = 1000) # location of your observation (observations x 2)
```

The method can run in the following way, where the number of signatures is set to 2 in this case:
```
res = nnmf(data = data, noSignatures = 2, location = location)

res$weights # The weights for each signature
res$signatures # The two signatures
```
### Multiple slices
If you have multiple slices, then you can simply set the id of the slice as the varaible batch. It has to be a vector of the same length as your observations.
```
slice_id = rep(c(1:10), each = 100) 
res = nnmf(data = data, noSignatures = 2, location = location, batch = slice_id)

res$weights # The weights for each signature
res$signatures # The two signatures

```
### Slices with more than 20.000 observations
If you have a very large dataset (more than 20.000 observations) it is recommended to split your dataset into batches. Use the groupondist function to do this. For large datasets you would set size between 10.000 to 20.000, where I here just set 100 as an example. 
```
batch_id = groupondist(location,size = 100)
```
Now you can run your model with the specified batches and you can start by using the predefined lengthscale, which you can later try to change after assessing your results.
```
res = nnmf(data = data, noSignatures = 2, location = location, batch = batch_id)

res$weights # The weights for each signature
res$signatures # The two signatures

```

### Find the top weighted features/genes in your data 
To find the top weighted genes of your signatures you can use the function topfeatures(). Below the top 10 genes are found from the data 
```
genes = paste0('gene',1:100)
top_name = topfeatures(signatures = res$signatures,feature_names = genes, ntop = 10)

```
## Citing the work 

For more details on the methods and examples go to the [Manuscript](https://www.biorxiv.org/content/10.1101/2025.04.26.650724v1).

