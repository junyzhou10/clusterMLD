#' @title Plot estimated mean trajectories for each detect cluster
#' @describeIn Based on output from LongDataCluster, yield the corresponding mean curves
#' @param Cluster.object The object return from function LongDataCluster
#' @param No.Cluster User specified number of clusters
#' @param CH If TRUE, iuse optimal number of clusters determined by CH index to cut dendrogram. The default is to use Gapb
#' @param add.sample The number of samples plot addtionally. Default is 0 means only plot mean trajectories
#' @param trsp Transparent coefficient when plotting samples
#' @param seed Seed to draw samples from the data
#' @param ... Additional arguments for plot
#' @import dplyr grDevices graphics
#' @examples \dontrun{
#' require(doMC)
#' registerDoMC(cores = 6)
#' output = LongDataCluster(Longdat2$Dat$obs,
#' Longdat2$Dat[,paste("y", seq(5), sep = "_")],
#' Longdat2$Dat$id, parallel = T, dropout = 15, part.size = 200)
#' MeanPlot(output, add.sample = 5)
#' }
#' @export


MeanPlot <- function(Cluster.object, No.Cluster = NULL, CH = FALSE, add.sample = 0, trsp = 0.5, seed = 123, ...) {
  args <- list(...)
  # inputs
  if (is.null(No.Cluster)) {
    NoCl  = ifelse(CH, Cluster.object$No.CH, Cluster.object$No.Gapb)
  } else {
    NoCl  = No.Cluster
  }

  # determine the number of columns according to the number of outcomes
  no.Y = dim(data.matrix(Cluster.object$calls$Y))[2]
  mod3 = no.Y %% 3
  mod4 = no.Y %% 4
  if (mod3==0) {
    ncolumn = 3
  } else if (mod4==0) {
    ncolumn = 4
  } else {
    ncolumn = ifelse(mod3<=mod4, 4, 3)
  }
  if (no.Y <= 4) {
    layout( matrix(seq(no.Y), nrow = 1, byrow = T) )
  } else {
    layout( matrix(seq(ceiling(no.Y/ncolumn)*ncolumn), nrow = ceiling(no.Y/ncolumn), byrow = T) )
  }


  # plot
  set.seed(seed)
  no.obs = dim(Cluster.object$calls$x.bs)[1]
  cluster.ID = Cluster.object$out.ID[[NoCl]]
  x.bs = Cluster.object$calls$x.bs
  id.seq  = Cluster.object$calls$id
  Y.dat = Cluster.object$calls$Y.dat
  id.sample = sample(no.obs, min(no.obs, 1000))
  x.sample  = Cluster.object$calls$x[id.sample]
  x.ord     = order(x.sample); x.sample = x.sample[x.ord]
  x.fit = Cluster.object$calls$x.bs[id.sample, ]; x.fit = x.fit[x.ord,]
  y.fitted = list()
  for (cc in seq(NoCl)) {
    x0  = x.bs[which(id.seq %in% cluster.ID[[cc]]), ]
    y0  = Y.dat[which(id.seq %in% cluster.ID[[cc]]), ]
    y.fitted = c(y.fitted, list(x.fit %*% ginv(t(x0) %*% x0) %*% t(x0) %*% t(t(y0))))
  }

  # prepare color for mean trajectories and samples
  col.list = rainbow(NoCl, start = 0.01, end = 0.8)
  col.list2 = apply(col2rgb(col.list), 2, function(x) rgb(red=x[1]/255,green=x[2]/255,blue=x[3]/255,alpha=trsp))

  for (j in seq(no.Y)) {
    tmp = data.frame(x=1,y=1); colnames(tmp)<-c("Observation", colnames(Y.dat)[j])
    plot(tmp, type="n", xlim=range(x.sample), ylim=range(y.fitted), ...)
    for (cc in seq(NoCl)) {
      lines(x.sample, y.fitted[[cc]][, j], col = col.list[cc], type = "l", lwd = 2)
    }

    if (add.sample > 0) {
      # randomly pick n=add.sample ids from each cluster
      id.sel = lapply(cluster.ID, function(x) sample(x, min(add.sample, length(x))))
      for (id in unlist(id.sel)) {
        lines(Cluster.object$calls$x[id.seq==id], Y.dat[id.seq==id, j], type = "b", col = col.list2[sapply(cluster.ID, function(x) id %in% x)], lwd = 1)
      }
    }
  }

}
