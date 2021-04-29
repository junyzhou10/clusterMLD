#' @title Plot estimated mean trajectories for each detected cluster
#' @description Based on output from LongDataCluster, yield the corresponding mean curves
#' @param Cluster.object The object return from function LongDataCluster
#' @param No.Cluster User specified number of clusters
#' @param CH If TRUE, use optimal number of clusters determined by CH index to cut dendrogram. The default is to use Gap.b
#' @return A figure yield by plotly, which mean trajectories in different colors.
#' @import dplyr ggplot2
#' @importFrom stats reshape
#' @importFrom plotly highlight ggplotly highlight_key
#' @examples \dontrun{
#' require(doMC)
#' registerDoMC(cores = 6)
#' output = LongDataCluster(Longdat2$Dat$obs,
#' Longdat2$Dat[,paste("y", seq(5), sep = "_")],
#' Longdat2$Dat$id, parallel = T, dropout = 15, part.size = 200)
#' MeanPlot(output)
#' }
#' @export


MeanPlot <- function(Cluster.object, No.Cluster = NULL, CH = FALSE) {
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
    Nrow = 1
  } else {
    Nrow = ceiling(no.Y/ncolumn)
  }


  # plot
  no.obs = dim(Cluster.object$calls$x.bs)[1]
  cluster.ID = Cluster.object$Cluster.Lists[[NoCl]]
  x.bs = Cluster.object$calls$x.bs
  id.seq  = Cluster.object$calls$id
  Y.dat = Cluster.object$calls$Y.dat
  id.sample = sample(no.obs, min(no.obs, 1000))
  x.sample  = Cluster.object$calls$x[id.sample];
  x.ord     = order(x.sample); x.sample = x.sample[x.ord]
  x.fit = Cluster.object$calls$x.bs[id.sample, ]; x.fit = x.fit[x.ord,]
  plot_dat = NULL
  for (cc in seq(NoCl)) {
    x0  = x.bs[which(id.seq %in% cluster.ID[[cc]]), ]
    y0  = Y.dat[which(id.seq %in% cluster.ID[[cc]]), ]
    plot_dat = rbind(plot_dat, data.frame(x.fit %*% ginv(t(x0) %*% x0) %*% t(x0) %*% t(t(y0)), x = x.sample, grp = paste("Cluster",cc)))
  }
  fitted = reshape(plot_dat, direction = "long", varying = colnames(Y.dat), v.names = "Y", timevar = "ylabel", times = colnames(Y.dat))

  p <- ggplot(highlight_key(fitted, key=~grp), aes(x, Y)) + geom_line(aes(color = grp)) +
    facet_wrap(~ylabel, nrow = Nrow) + theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Observations", y = "Outcomes")
  gg <- highlight(ggplotly(p, tooltip = c("grp")), "plotly_hover", 'plotly_doubleclick')
  return(gg)
}
