#' @title Clustering longitudinal data (Main Function)
#' @author Junyi Zhou <\email{junyzhou@iu.edu}>
#' @description Clustering longitudinal data, expecially tailored for those observational longitudinal data with sparse and irregular observations. The output could be a vector, i.e., at each occasion, more than one measure are observed for each subject
#' @param x A vector in long format, occassions or time of observation times
#' @param Y A matrix, if multiple outcomes exist; or a (column) vector, for single outcome case
#' @param id A vector with same length of x, represents the corresponding subject id of each observation
#' @param DistMetric \code{c("W", "UnW")} for Weighted or Unweighted distance metric. Default is \code{Weighted}. See details.
#' @param functional A string from \code{c("bs", "ns")}, indicating b-splines or natural splines
#' @param preprocess boolean, whether data pre-processing procedure should be applied. Default is \code{TRUE} to handle 
#'                   subjects with observations less than number of parameters. If set to \code{FALSE}, those subjects will be excluded
#' @param weight.func A string from \code{c("standardize", "softmax", "equal")}, a method to handle weights across multiple outcomes. 
#'                    Default is \code{"standardize"}, but \code{"softmax"} is recommended for cases with a lot noise (indistinguishable) outcomes.
#'                    Select \code{"equal"} then all outcomes are forced to be treated equally.
#' @param parallel If \code{TRUE}, use parallel foreach to cluster each ubgroup of original data. 
#'                 Must register parallel beforehand, such as doMC or others
#' @param stop An integer, only required in parallel computing. It indicates when to stop hierarchical 
#'             clustering process in each slave/server. See details
#' @param part.size In parallel computing, the (rough) number of subjects in each random partitions. See more in details
#' @param ... Additional arguments for the functional bases. The default is cubic B-spline with 3 interval knots
#' @details There are two distance metrics available, one is so-called unweighted \code{UnW}, i.e. \cr
#'          \eqn{Dist = SSR(both) - SSR(cluster 1) - SSR(cluster 2)}. 
#'          \cr The other is weighted \code{W}, \cr
#'          \eqn{Dist = {SSR(both) - SSR(cluster 1) - SSR(cluster 2)}/p / {(SSR(cluster 1) + SSR(cluster 2))/(n_1+n_2-p)}},
#'          \cr where \eqn{p} is the number of basis function, and \eqn{n_1} and \eqn{n_2} are the sample size for cluster 1
#'          and 2, respectively. Two metric can yield slightly different results. 
#'          \cr \cr
#'          For relatively large sample size, in terms of the number of subjects not observations, we suggest to apply parallel computing to save time dramatically.
#'          By specifying \code{part.size} and \code{stop}, the algorithm actually split data into multiple random partitions with size roughly 
#'          equal to part.size, and then apply the hierarchical algorithm in a parallel fashion on each partition until the number of clusters goes down to \code{stop}.
#'          Then combine the output clusters together. If the remaining number of clusters is still larger than part.size, multiple rounds of this process will be implemented.
#' @return A list object containing the hierarchical clustering results, and some ancillary outputs for parallel computing.
#' \item{Cluster.res}{A list with length equal to the number of clusters determined by \eqn{Gap_b} index. Each list element consists of the corresponding subjects' id. In particular, it is equivalent to \code{out.ID[[No.Gapb]]}}
#' \item{Cluster.Lists}{List of hierarchical results where the length should be the number of subjects when \code{parallel = FALSE}. This is the main output. The rest elements in returns are used to facilitate other functions}
#' \item{No.Gapb}{Number of clusters determined by \eqn{Gap_b} statistic}
#' \item{No.CH}{Number of clusters determined by CH index}
#' @import splines foreach
#' @importFrom MASS ginv
#' @examples
#' output = LongDataCluster(Longdat$Dat$obs,
#' Longdat$Dat[,paste("y", seq(5), sep = "_")], Longdat$Dat$id)
#' # parallel version
#' \dontrun{
#' require(doMC)
#' registerDoMC(cores = 6)
#' output = LongDataCluster(Longdat2$Dat$obs,
#' Longdat2$Dat[,paste("y", seq(5), sep = "_")],
#' Longdat2$Dat$id, parallel = T, stop = 15, part.size = 200)
#' }
#' @export
#' @seealso \code{\link{MeanPlot}}, \code{\link{DendroPlot}}

LongDataCluster <- function(x, Y, id, 
                            DistMetric = "W",
                            functional = "bs", 
                            preprocess = TRUE, 
                            weight.func = "standardize", 
                            parallel = FALSE, stop = 20, part.size = 300, ...){
  id.list = unique(id)
  y.dim = ncol(data.matrix(Y))
  if (parallel) {
    splits = floor(length(unique(id))/part.size)
    # randomly split data into subgroups
    id.split = sample(rep(seq(splits), length.out = length(id.list)))
    pure.leaves = foreach(ii = seq(splits), .combine = c) %dopar% {
      grpID = id.list[id.split==ii]
      if (DistMetric == "W") {
        output = LongDataClusterW(x=x[id %in% grpID], Y=data.matrix(Y)[id %in% grpID,], id=id[id %in% grpID], functional = functional, preprocess = preprocess, weight.func = weight.func, parallel = TRUE, stop = stop, ...)
      } else if (DistMetric == "UnW") {
        output = LongDataClusterUnW(x=x[id %in% grpID], Y=data.matrix(Y)[id %in% grpID,], id=id[id %in% grpID], functional = functional, preprocess = preprocess, weight.func = weight.func, parallel = TRUE, stop = stop, ...)
      }
      
      return(output$pure.leaf)
    }

    # iteratively apply parallel computing idea
    while (TRUE) {
      len = length(pure.leaves)
      if (len <= part.size*1.4) {
        break
      }
      # recalculate how many partitions required
      splits = ceiling(len/part.size)
      id.split = sample(rep(seq(splits), length.out = len))
      pure.leaves = foreach(ii = seq(splits), .combine = c) %dopar% {
        output = FinalCluster(pure.leaf = pure.leaves[id.split == ii], stop = stop, weight.func = weight.func, DistMetric = DistMetric, y.dim = y.dim)
        return(output$pure.leaf)
      }
    }
    res = FinalCluster(pure.leaf = pure.leaves, weight.func = weight.func, DistMetric = DistMetric, y.dim = y.dim) # stop = 1, combine to 1 group
  } else {
    if (DistMetric == "W") {
      res = LongDataClusterW(x=x, Y=Y, id=id, functional = functional, preprocess = preprocess, weight.func = weight.func, ...)
    } else if (DistMetric == "UnW") {
      res = LongDataClusterUnW(x=x, Y=Y, id=id, functional = functional, preprocess = preprocess, weight.func = weight.func, ...)
    }
  }

  # final output
  args <- list(...)
  if (!"df" %in% names(args)) {
    if (functional == "bs") {
      x.bs = bs( x, df = 7, intercept = T, ...)
    } else if (functional == "ns") {
      x.bs = ns( x, df = 7, intercept = T, ...)
    }

  } else {
    if (functional == "bs") {
      x.bs = bs( x, intercept = T, ...)
    } else if (functional == "ns") {
      x.bs = ns( x, intercept = T, ...)
    }
  }

  Gap_b     = res$Gap.dist - seq(1, length(res$Gap.dist))*res$p.var*res$e.sigma
  # Gap_b     = res$Gap.dist - seq(1, length(res$Gap.dist))*(length(x) - 2*res$p.var)/(length(x) - 2*res$p.var - 2) # when using Chow-test as distance
  N.cluster = ifelse(sum(diff(Gap_b)<0)==0, 1, which(diff(Gap_b)<0)[1])
  CH.index  = res$B.dist/(seq(1,length(res$B.dist))-1)/res$W.dist*(sum(res$obs.no) - res$p.var*seq(1,length(res$B.dist))); CH.index[1] = 0
  return(list(Cluster.res = lapply(res$out.ID[[N.cluster]], sort),
              Cluster.Lists = res$out.ID,
              No.Gapb= N.cluster,
              No.CH  = which.max(CH.index),
              Gap_b  = Gap_b,
              CH.index = CH.index,
              weight = res$weight,
              Name.L = res$Name.L,
              Addres = res$Addres,
              calls  = list(x.bs = x.bs,
                            x = x,
                            Y.dat = data.matrix(Y),
                            id= id,
                            functional = functional,
                            preprocess = preprocess,
                            weight.func= weight.func,
                            parallel   = parallel,
                            stop       = stop,
                            part.size  = part.size)))
}

