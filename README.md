# ClusterCall

ClusterCall is an R package for making tetraploid genotype calls from SNP array data. This repository was created to facilitate the instalation of ClusterCall during the EiB workshop on Polyploid Data management and analysis CIP, Peru May 8-10. This package was developed by Jeff Endelman Group . The original program source can be accessed [here](https://potatobreeding.cals.wisc.edu/wp-content/uploads/sites/161/2017/08/ClusterCall_Download.zip) and the raleted publication can be found in [Schmitz Carley et al. 2017](https://doi.org/10.1007/s00122-016-2845-5).

ClusterCall is not available in CRAN, but you can install it from Git Hub. Within R, you need to install and load the package `devtools`:

```R
install.packages("devtools")
```
To install ClusterCall from Git Hub use

```R
devtools::install_github("mmollina/ClusterCall")
```