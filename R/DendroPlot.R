#' @title Plot dendrogram
#' @describeIn Based on output from LongDataCluster, yield the corresponding dendropgram
#' @param Cluster.object The object return from function LongDataCluster
#' @param No.Cluster User specified number of clusters
#' @param CH If TRUE, iuse optimal number of clusters determined by CH index to cut dendrogram. The default is to use Gapb
#' @param plain If TRUE, the dendrogram only, without coloring the clustering results
#' @param ... Additional arguments for plot
#' @examples
#' output = LongDataCluster(Longdat$Dat$obs,
#' Longdat$Dat[,paste("y", seq(5), sep = "_")],
#' Longdat$Dat$id)
#' DendroPlot(output)
#'
#' @export


DendroPlot <- function(Cluster.object, No.Cluster = NULL, CH = FALSE, plain = FALSE, ...){
  args <- list(...)
  # inputs
  if (is.null(No.Cluster)) {
    NoCl  = ifelse(CH, Cluster.object$No.CH, Cluster.object$No.Gapb)
  } else {
    NoCl  = No.Cluster
  }

  Add.res = rev(Cluster.object$Addres)
  K       = length(Cluster.object$Name.L)
  Name.L  = Cluster.object$Name.L

  Level = K + 1 - NoCl
  Add.res  = log(Add.res) - min(log(Add.res)) + 0.1
  Grp.Size = seq(K)
  col.list = rainbow(NoCl, start = 0.01, end = 0.8)



  # Figure Frame
  plot(data.frame(Clusters=0,Height = 0), type = "n", ylim = c(0,max(Add.res)+0.1),xlim = c(0, K), axes = FALSE, ...)
  box()

  # Indicating CLUSTER by color
  ord.knot = as.numeric(unlist(strsplit(unlist(Name.L[[length(Name.L)]]),split=',')))
  pos.range = NULL
  for (i in 1:length(Name.L[[Level]])) {
    pos = which(ord.knot %in% as.numeric(unlist(strsplit(Name.L[[Level]][i], ","))))
    pos.range = rbind(pos.range, range(pos))
  }

  x.knot   = seq(1,K)[order(ord.knot)]
  y.knot   = array(0, K)
  a        = as.character(seq(1,K))


  for (k in 1:length(Add.res)) {
    b = Name.L[[k+1]]
    x = x.knot[!a %in% b]
    y = y.knot[!a %in% b]
    grpsize = Grp.Size[!a %in% b]
    if (k >= Level) {
      segments(x[1], y[1], x[1], Add.res[k], col= "grey45")
      segments(x[2], y[2], x[2], Add.res[k], col= "grey45")
      segments(x[1], Add.res[k], x[2], Add.res[k], col= "grey45")
    } else {
      if (plain){
        col.id = "grey45"
      } else {
        col.id = col.list[apply(pos.range, 1, function(rr) rr[1]<=x[1] & rr[2] >= x[1])]
      }

      segments(x[1], y[1], x[1], Add.res[k], col = col.id)
      segments(x[2], y[2], x[2], Add.res[k], col = col.id)
      segments(x[1], Add.res[k], x[2], Add.res[k], col = col.id)
    }

    x.knot = c(x.knot[a %in% b], mean(x))
    y.knot = c(y.knot[a %in% b], Add.res[k])
    Grp.Size = c(Grp.Size[a %in% b], sum(grpsize))
    a = b
  }

  if (!plain){
    abline(a=mean(Add.res[(Level-1):(Level)]), b = 0, col='gray29',lty = 2)
  }

}

