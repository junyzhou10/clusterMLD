#' @title Ancillary function for parallel computing (not for call)
#' @description Based on outputs from each server/slave, merge them and keep hierarchical algorithm until only one group left
#' @param pure.leaf Output from LongDataCluster.core in each server
#' @param stop When iteratively applying parallel idea, stop when a certain number of clusters are left
#' @param weight.func weight function
#' @param DistMetric Weighted or Unweighted distance metric
#' @return The same output as LongDataCluster.core

FinalCluster <- function(pure.leaf, stop = 1, weight.func, DistMetric) {
  # define weight function
  if (weight.func == "standardize") {
    w.func <- function(w) {w/sum(w)}
  } else if (weight.func == "softmax") {
    w.func <- function(w) {w = scale(w); exp(w)/sum(exp(w))}
  } else if (weight.func == "equal") {
    w.func <- function(w) {rep(1/y.dim, y.dim)}
  } else {
    warning("weight.func is not correctly specified!")
    return(NULL)
  }


  XY   = Reduce("+", lapply(pure.leaf, function(x) x$Xy))
  XX   = Reduce("+", lapply(pure.leaf, function(x) x$XX))
  Y2   = Reduce("+", lapply(pure.leaf, function(x) x$Y2))
  e.beta  = ginv(XX) %*% XY
  SSR.all = Y2 - diag(t(XY) %*% ginv(XX) %*% XY)
  obs.no = unlist(lapply(pure.leaf, function(x) x$N.in))
  p.var = dim(XX)[1]
  e.sigma = SSR.all/(sum(obs.no) - p.var)


  ## generate weight
  weight = SSR.all/(colSums(do.call(rbind, lapply(pure.leaf, function(x) x$SSR0))))-1;
  weight = w.func(weight)

  ## Generate all inter-cluster distances
  Dist.Gap = NULL # inter-cluster (or so-called between cluster distance)
  Dist.inter = NULL # between cluster distance for CH
  for (i in seq(length(pure.leaf))) {
    xy.s = XY - pure.leaf[[i]]$Xy
    xx.s = XX - pure.leaf[[i]]$XX
    SSR.s= Y2 - pure.leaf[[i]]$Y2 - diag(t(xy.s) %*% ginv(xx.s) %*% xy.s)
    Dist.Gap = c(Dist.Gap, sum((SSR.all - SSR.s - pure.leaf[[i]]$SSR0)*weight))
    # Dist.Gap = c(Dist.Gap, sum((SSR.all - SSR.s - pure.leaf[[i]]$SSR0)/(SSR.s + pure.leaf[[i]]$SSR0)*(sum(obs.no) - 2*p.var)/p.var*weight))
    
    e.beta.s = e.beta - ginv(pure.leaf[[i]]$XX) %*% pure.leaf[[i]]$Xy 
    Dist.inter = c(Dist.inter, sum(t(e.beta.s) %*% pure.leaf[[i]]$XX %*% e.beta.s *weight))
  }

  ## generate Dist.tab
  n.leaf   = length(pure.leaf)
  Dist.tab = matrix(Inf, n.leaf, n.leaf)

  # An original distance matrix to record all pairewise distance between pure subgroups
  for ( i in seq(n.leaf-1) ) {
    for ( j in (i+1):n.leaf ) {
      Y2.merge = pure.leaf[[i]]$Y2 + pure.leaf[[j]]$Y2
      Xy.merge = pure.leaf[[i]]$Xy + pure.leaf[[j]]$Xy
      XX.merge = pure.leaf[[i]]$XX + pure.leaf[[j]]$XX
      A.merge  = ginv(XX.merge)
      SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
      if (DistMetric == "UnW") {
        Dist.tab[i,j] = sum((SSR.merge - pure.leaf[[i]]$SSR0 - pure.leaf[[j]]$SSR0)*weight) # if knots changed during pure group generation, this should be SSR.merge/(df-dfi-dfj)
      } else if (DistMetric == "W") {
        Dist.tab[i,j] = sum(((SSR.merge - pure.leaf[[i]]$SSR0 - pure.leaf[[j]]$SSR0)/(pure.leaf[[i]]$SSR0 + pure.leaf[[j]]$SSR0)*(pure.leaf[[i]]$N.in + pure.leaf[[j]]$N.in - 2*p.var))*weight)
      }
    }
  }

  ##============= Based on pure.leaf and Dist.tab, start hierachical merging ==============##
  out.ID  = list(lapply(pure.leaf, function(x) x$id.in))
  B.dist  = sum(Dist.inter)
  Name.L   = list(as.character(seq(length(pure.leaf))))
  Addres   = numeric(0)
  W.dist   = sum(colSums(do.call(rbind, lapply(pure.leaf, function(x) x$SSR0)))*weight)
  Gap.dist = sum(Dist.Gap)

  while (TRUE) {
    ij = which(Dist.tab==min(Dist.tab), arr.ind = T)
    i  = ij[1]
    j  = ij[2]

    # 1. First update pure.leaf by removing ij and adding one combined subgroup
    Y2.merge = pure.leaf[[i]]$Y2 + pure.leaf[[j]]$Y2
    Xy.merge = pure.leaf[[i]]$Xy + pure.leaf[[j]]$Xy
    XX.merge = pure.leaf[[i]]$XX + pure.leaf[[j]]$XX
    # A.merge  = pure.leaf[[i]]$A.mat - pure.leaf[[i]]$A.mat %*% ginv(pure.leaf[[i]]$A.mat + pure.leaf[[j]]$A.mat) %*% pure.leaf[[i]]$A.mat
    A.merge  = ginv(XX.merge)
    SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
    id.merge = c(pure.leaf[[i]]$id.in, pure.leaf[[j]]$id.in)
    N.merge  = pure.leaf[[i]]$N.in + pure.leaf[[j]]$N.in
    leaf.merge = list(id.in = id.merge,
                      A.mat = A.merge,
                      XX    = XX.merge,
                      Xy    = Xy.merge,
                      Y2    = Y2.merge,
                      N.in  = N.merge,
                      SSR0  = SSR.merge)

    # farthest distance / second smallest distance

    # B.mean = c(mean(tmp.D), B.mean)
    # B.sd   = c(sd(tmp.D), B.sd)
    pq = which(Dist.tab==max(Dist.tab[is.finite(Dist.tab)]), arr.ind = T)
    ii  = pq[1]
    jj  = pq[2]

    # update leaf info, SSR, and ID
    pure.leaf = c(pure.leaf[-ij], list(leaf.merge))
    W.dist    = c(sum(colSums(do.call(rbind, lapply(pure.leaf, function(x) x$SSR0)))*weight), W.dist)
    out.ID    = c(list(c(out.ID[[1]][-ij], list(unlist(out.ID[[1]][ij])))), out.ID)
    Name.L    = c(Name.L, list( c(Name.L[[length(Name.L)]][-ij],  paste0(c(Name.L[[length(Name.L)]][i], Name.L[[length(Name.L)]][j]), collapse = ","))  ))
    Addres    = c(min(Dist.tab)/p.var, Addres)

    # update Dist.inter
    xy.s = XY - leaf.merge$Xy
    xx.s = XX - leaf.merge$XX
    SSR.s= Y2 - leaf.merge$Y2 - diag(t(xy.s) %*% ginv(xx.s) %*% xy.s)
    Dist.Gap = c(Dist.Gap[-ij], sum(( SSR.all - SSR.s - leaf.merge$SSR0)*weight))
    
    e.beta.s = e.beta - ginv(leaf.merge$XX) %*% leaf.merge$Xy 
    Dist.inter = c(Dist.inter[-ij], sum(t(e.beta.s) %*% leaf.merge$XX %*% e.beta.s *weight))
    # Dist.inter = c(Dist.inter[-ij], sum(( SSR.all - SSR.s - leaf.merge$SSR0)/(SSR.s + leaf.merge$SSR0)*(sum(obs.no) - 2*p.var)/p.var*weight))
    B.dist  = c(sum(Dist.inter), B.dist)
    Gap.dist = c(sum(Dist.Gap), Gap.dist)
    
    
    n.tmp = length(pure.leaf)
    # 2. Update distance matrix by removing ij, and add one row/column with distance between new subgroup to all rest subgroups
    Dist.new = NULL
    Dist.tab = data.matrix(Dist.tab[-ij, -ij]);

    if (n.tmp == 1) { # no need for further combine
      Dist.tab = leaf.merge$SSR0
      break
    } else {
      for ( k in seq(dim(Dist.tab)[1]) ) {
        Y2.merge = pure.leaf[[k]]$Y2 + leaf.merge$Y2
        Xy.merge = pure.leaf[[k]]$Xy + leaf.merge$Xy
        XX.merge = pure.leaf[[k]]$XX + leaf.merge$XX
        # A.merge  = pure.leaf[[k]]$A.mat - pure.leaf[[k]]$A.mat %*% ginv(pure.leaf[[k]]$A.mat + leaf.merge$A.mat) %*% pure.leaf[[k]]$A.mat
        A.merge  = ginv(XX.merge)
        SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
        if (DistMetric == "UnW") {
          Dist.new = c(Dist.new, sum(weight*(SSR.merge - pure.leaf[[k]]$SSR0 - leaf.merge$SSR0)))
        } else if (DistMetric == "W") {
          Dist.new = c(Dist.new, sum(weight*((SSR.merge - pure.leaf[[k]]$SSR0 - leaf.merge$SSR0)/(pure.leaf[[k]]$SSR0 + leaf.merge$SSR0)*(pure.leaf[[k]]$N.in + leaf.merge$N.in - 2*p.var))))
        }
      }
      Dist.tab = cbind(Dist.tab, Dist.new)
      Dist.tab = rbind(Dist.tab, rep(Inf, dim(Dist.tab)[2]))
    }

    if (n.tmp <= stop) {
      break
    }
  }

  return(list(out.ID = out.ID,
              pure.leaf = pure.leaf,
              Dist.tab = Dist.tab,
              B.dist = B.dist,
              W.dist = W.dist,
              Gap.dist = Gap.dist,
              obs.no = obs.no,
              weight = weight,
              Name.L = Name.L,
              Addres = Addres,
              p.var = p.var,
              e.sigma= sum(weight*e.sigma)))
}


