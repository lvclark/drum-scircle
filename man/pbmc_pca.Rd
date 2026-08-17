\name{pbmc_pca}
\alias{pbmc_pca}
\docType{data}
\title{
Principal Components Analysis and cluster assignment from 3K PBMCs
}
\description{
This testing dataset is generated using the Seurat tutorial from the 3K PBMC
dataset provided by 10X genomics. It includes a matrix of 50 principal components
across 2638 cells, as well as a factor indicating cluster assignments to those
cells.
}
\usage{data("pbmc_pca")}
\format{
  The format is:
 num [1:2638, 1:50] -4.73 -0.518 -3.189 12.793 -3.129 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:2638] "AAACATACAACCAC-1" "AAACATTGAGCTAC-1" "AAACATTGATCAGC-1" "AAACCGTGCTTCCG-1" ...
  ..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...

  and:
 Factor w/ 9 levels "0","1","2","3",..: 3 4 3 2 7 3 5 5 5 6 ...
 - attr(*, "names")= chr [1:2638] "AAACATACAACCAC-1" "AAACATTGAGCTAC-1" "AAACATTGATCAGC-1" "AAACCGTGCTTCCG-1" ...
}
\source{
https://satijalab.org/seurat/articles/pbmc3k_tutorial
}
\examples{
data(pbmc_pca)

table(pbmc_clust)

pbmc_pca[1:10,1:10]
}
\keyword{datasets}
