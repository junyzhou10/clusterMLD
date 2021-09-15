#' @title Clustering longitudinal data (not for call)
#' @description Clustering longitudinal data, corresponding to weighted DistMetric argument,, expecially tailored for those observational longitudinal data with sparse and irregular observations. The output could be a vector, i.e., at each occasion, more than one measure are observed for each subject.
#' @param x A vector in long format, occassions or time of observation times.
#' @param Y A matrix, if multiple outcomes exist; or a (column) vector, for single outcome case
#' @param id A vector with same length of x, represents the corresponding subject id of each observation
#' @param functional A string from c("bs", "ns"), indicating b-splines and natural splines
#' @param preprocess boolean, whether data pre-processing procedure should be applied. Default is TRUE to handle subjects with observations less than number of parameters. If set to FALSE, those subjects will be excluded.
#' @param weight.func A string from c("standardize", "softmax"), a method to handle weights across multiple outcomes. Default is "standardize", but "softmax" is recommended for cases with a lot noise (indistinguishable) outcomes
#' @param parallel boolean, indicating whether parallel computing should be implemented. Default is FALSE
#' @param stop An integer, only required in parallel computing. It indicates when to stop hierarchical clustering process in each slave/serve
#' @param ... Additional arguments for the functional bases. The default is cubic B-spline with 3 interval knots
#' @return A list object containing the hierarchical clustering results, and some ancillary outputs for parallel computing. The optimal number of clusters will not be determined by this function.

LongDataClusterW <- function(x, Y, id, functional = "bs", preprocess = TRUE, weight.func = "standardize", parallel = FALSE, stop = 20, ...) {
  # define weight function
  if (weight.func == "standardize") {
    w.func <- function(w) {w/sum(w)}
  } else if (weight.func == "softmax") {
    w.func <- function(w) {w = scale(w); exp(w)/sum(exp(w))}
  } else {
    warning("weight.func is not correctly specified!")
    return(NULL)
  }
  
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
  
  x.list = split(as.data.frame(x.bs),  id)
  x.list = lapply(x.list, data.matrix)
  Y.dat  = apply(data.matrix(Y), 2, scale) # in case input Y is a vector
  y.list = split(as.data.frame(Y.dat),  id)
  y.list = lapply(y.list, data.matrix)
  
  # some global variables
  p.var   = ncol(x.bs) # number of parameters
  y.dim   = ncol(Y.dat) # number of outcomes
  id.list = id.wait = unique(id) # individual id
  id.seq  = id # total id sequence, used a lot in filter data
  obs.no  = unlist(lapply(x.list, function(x) nrow(x))) # the # observations for each subject
  
  x.wait  = x.bs[which(id.seq %in% id.wait), ]
  y.wait  = data.matrix(Y.dat[which(id.seq %in% id.wait),])
  XX.wait = XX.remain = t(x.wait) %*% x.wait
  Xy.wait = Xy.remain = t(x.wait) %*% t(t(y.wait))
  Y2.wait = Y2.remain = colSums(y.wait^2)
  A.wait  = ginv(XX.wait)
  e.beta  = A.wait %*% Xy.wait
  SSR.all = Y2.wait - diag(t(Xy.wait) %*% e.beta)
  e.sigma = SSR.all/(nrow(x.wait) - p.var)
  
  ## Generate profile for each subject
  pure.leaf = NULL
  for (i in seq(length(x.list))) {
    x.in    = x.list[[i]]
    y.in    = y.list[[i]]
    N.in    = nrow(x.in)
    XX      = t(x.in) %*% x.in
    A.mat   = ginv(XX)
    Xy      = t(x.in) %*% t(t(y.in))
    Y2      = colSums(y.in^2)
    SSR0    = diag(Y2 - t(Xy) %*% A.mat %*% Xy)
    SSR.r   = diag(Y2.wait - Y2 - t(Xy.wait-Xy) %*% solve(XX.wait - XX) %*% (Xy.wait-Xy))
    # SSR.dist= SSR.all - SSR0 - SSR.r
    SSR.dist= (SSR.all - SSR0 - SSR.r)/(SSR0 + SSR.r)*(sum(obs.no) - 2*p.var) # p could be neglected
    leaf.new = list(id.in = id.list[i], A.mat = A.mat, XX = XX, Xy = Xy, Y2 = Y2, N.in = N.in, SSR0 = SSR0, SSR.dist = SSR.dist)
    
    pure.leaf = c(pure.leaf, list(leaf.new))
  }
  pure.leaf0 = pure.leaf
  
  
  if (preprocess) {
    ## Check if all subjects have enough observations. If not, implement data pre-preprocess
    if (min(obs.no) <= p.var) { # some of them do not have enough observations, need further combination of subjects
      if (max(obs.no) <= p.var) { # equally assign weight
        weight = array(1/y.dim, y.dim)
      } else {
        id.enough = id.list[which(obs.no > p.var)]
        x0  = x.bs[which(id.seq %in% id.enough), ]
        y0  = data.matrix(Y.dat[which(id.seq %in% id.enough),])
        XX0 = t(x0) %*% x0
        Xy0 = t(x0) %*% t(t(y0))
        Y20 = colSums(y0^2)
        ssr0 = Y20 - diag(t(Xy0) %*% ginv(XX0) %*% Xy0)
        
        weight = ssr0/(colSums(do.call(rbind, lapply(pure.leaf[which(obs.no > p.var)], function(x) x$SSR0))))-1;
        weight = w.func(weight)
      }
      ##================== based on weight, start to combine those subject with insufficient observations ==============##
      # the combination of individual is determined by the closest MSR
      # we first generate the MSR table for all pairs, if after combination still could not calculate SSR, leave Inf
      n.leaf   = length(pure.leaf)
      Dist.tab = matrix(Inf, n.leaf, n.leaf)
      
      # An original distance matrix to record all pairewise distance between pure subgroups
      for ( i in seq(n.leaf-1) ) {
        for ( j in (i+1):n.leaf ) {
          N.in     = pure.leaf[[i]]$N.in + pure.leaf[[j]]$N.in
          if (N.in > p.var) {
            Y2.merge = pure.leaf[[i]]$Y2 + pure.leaf[[j]]$Y2
            Xy.merge = pure.leaf[[i]]$Xy + pure.leaf[[j]]$Xy
            XX.merge = pure.leaf[[i]]$XX + pure.leaf[[j]]$XX
            A.merge  = ginv(XX.merge)
            SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
            Dist.tab[i,j] = sum(SSR.merge*weight)/N.in
          }
        }
      }
      
      # then, combine the individuals with smallest MSR, and update the distance matrix
      while( TRUE ) {
        length.in = unlist(lapply(pure.leaf, function(x) x$N.in))
        if (sum(length.in <= p.var) < 1) {
          break
        }
        
        Dist.tmp = Dist.tab[which(length.in <= p.var), ]
        
        ij = which(Dist.tmp==min(Dist.tmp), arr.ind = T)
        if (length(ij) == 1) {
          ij = c(1, ij)
        }
        i  = ij[1]; i = which(length.in <= p.var)[i]
        j  = ij[2]
        
        
        # 1. First update pure.leaf by removing ij and adding one combined subgroup
        Y2.merge = pure.leaf[[i]]$Y2 + pure.leaf[[j]]$Y2
        Xy.merge = pure.leaf[[i]]$Xy + pure.leaf[[j]]$Xy
        XX.merge = pure.leaf[[i]]$XX + pure.leaf[[j]]$XX
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
        # update leaf info, SSR, and ID
        pure.leaf = c(pure.leaf[-c(i,j)], list(leaf.merge))
        # !!!!!!    update weight     !!!!!
        length.in = unlist(lapply(pure.leaf, function(x) x$N.in))
        id.enough    = unlist(lapply(pure.leaf[length.in > p.var], function(x) x$id.in))
        x0  = x.bs[which(id.seq %in% id.enough), ]
        y0  = data.matrix(Y.dat[which(id.seq %in% id.enough),])
        XX0 = t(x0) %*% x0
        Xy0 = t(x0) %*% t(t(y0))
        Y20 = colSums(y0^2)
        ssr0 = Y20 - diag(t(Xy0) %*% ginv(XX0) %*% Xy0)
        
        
        # notice after the first merge, weight will be 0 if using previous calculation, so still treat it as equal weight
        if (length(pure.leaf) < (n.leaf-1) ) {
          weight    = ssr0/(colSums(do.call(rbind, lapply(pure.leaf[which(length.in > p.var)], function(x) x$SSR0))))-1;
          # weight[3:5] = weight[3:5]/3
          weight = w.func(weight)
        }
        
        
        # 2. Update distance matrix by removing ij, and add one row/column with distance between new subgroup to all rest subgroups
        Dist.new = NULL
        Dist.tab = data.matrix(Dist.tab[-c(i,j), -c(i,j)]);
        
        for ( k in seq(nrow(Dist.tab)) ) {
          Y2.merge = pure.leaf[[k]]$Y2 + leaf.merge$Y2
          Xy.merge = pure.leaf[[k]]$Xy + leaf.merge$Xy
          XX.merge = pure.leaf[[k]]$XX + leaf.merge$XX
          N.merge  = pure.leaf[[k]]$N.in + leaf.merge$N.in
          A.merge  = ginv(XX.merge)
          SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
          Dist.new = c(Dist.new, sum(weight*SSR.merge)/N.merge)
        }
        Dist.tab = cbind(Dist.tab, Dist.new)
        Dist.tab = rbind(Dist.tab, rep(Inf, ncol(Dist.tab)))
      }
    } else {
      weight = SSR.all/(colSums(do.call(rbind, lapply(pure.leaf, function(x) x$SSR0))))-1
      weight = w.func(weight)
      n.leaf   = length(pure.leaf)
      Dist.tab = matrix(Inf, n.leaf, n.leaf)
      
      # An original distance matrix to record all pairewise distance between pure subgroups
      for ( i in seq(n.leaf-1) ) {
        for ( j in (i+1):n.leaf ) {
          Y2.merge = pure.leaf[[i]]$Y2 + pure.leaf[[j]]$Y2
          Xy.merge = pure.leaf[[i]]$Xy + pure.leaf[[j]]$Xy
          XX.merge = pure.leaf[[i]]$XX + pure.leaf[[j]]$XX
          # A.merge  = pure.leaf[[i]]$A.mat - pure.leaf[[i]]$A.mat %*% ginv(pure.leaf[[i]]$A.mat + pure.leaf[[j]]$A.mat) %*% pure.leaf[[i]]$A.mat
          A.merge  = ginv(XX.merge)
          SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
          # Dist.tab[i,j] = sum((SSR.merge - pure.leaf[[i]]$SSR0 - pure.leaf[[j]]$SSR0)*weight) # if knots changed during pure group generation, this should be SSR.merge/(df-dfi-dfj)
          Dist.tab[i,j] = sum(((SSR.merge - pure.leaf[[i]]$SSR0 - pure.leaf[[j]]$SSR0)/(pure.leaf[[i]]$SSR0 + pure.leaf[[j]]$SSR0)*(pure.leaf[[i]]$N.in + pure.leaf[[j]]$N.in - 2*p.var))*weight)
        }
      }
    }
  } else {
    id.enough = id.list[which(obs.no > p.var)]
    x0  = x.bs[which(id.seq %in% id.enough), ]
    y0  = data.matrix(Y.dat[which(id.seq %in% id.enough),])
    XX0 = t(x0) %*% x0
    Xy0 = t(x0) %*% t(t(y0))
    Y20 = colSums(y0^2)
    ssr0 = Y20 - diag(t(Xy0) %*% ginv(XX0) %*% Xy0)
    
    weight = ssr0/(colSums(do.call(rbind, lapply(pure.leaf[which(obs.no > p.var)], function(x) x$SSR0))))-1;
    weight = w.func(weight)
    pure.leaf = pure.leaf[which(obs.no > p.var)]
    n.leaf   = length(pure.leaf)
    Dist.tab = matrix(Inf, n.leaf, n.leaf)
    
    # An original distance matrix to record all pairewise distance between pure subgroups
    for ( i in seq(n.leaf-1) ) {
      for ( j in (i+1):n.leaf ) {
        Y2.merge = pure.leaf[[i]]$Y2 + pure.leaf[[j]]$Y2
        Xy.merge = pure.leaf[[i]]$Xy + pure.leaf[[j]]$Xy
        XX.merge = pure.leaf[[i]]$XX + pure.leaf[[j]]$XX
        # A.merge  = pure.leaf[[i]]$A.mat - pure.leaf[[i]]$A.mat %*% ginv(pure.leaf[[i]]$A.mat + pure.leaf[[j]]$A.mat) %*% pure.leaf[[i]]$A.mat
        A.merge  = ginv(XX.merge)
        SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
        # Dist.tab[i,j] = sum((SSR.merge - pure.leaf[[i]]$SSR0 - pure.leaf[[j]]$SSR0)*weight) # if knots changed during pure group generation, this should be SSR.merge/(df-dfi-dfj)
        Dist.tab[i,j] = sum((SSR.merge - pure.leaf[[i]]$SSR0 - pure.leaf[[j]]$SSR0)/(pure.leaf[[i]]$SSR0 + pure.leaf[[j]]$SSR0)*(pure.leaf[[i]]$N.in + pure.leaf[[j]]$N.in - 2*p.var)*weight)
      }
    }
  }
  
  
  
  ## Generate all inter-cluster distances for Gap_b calculation
  Dist.Gap = NULL # inter-cluster (or so-called between cluster distance)
  Dist.inter = NULL # between cluster distance for CH
  for (i in seq(length(pure.leaf))) {
    xy.s = Xy.wait - pure.leaf[[i]]$Xy
    xx.s = XX.wait - pure.leaf[[i]]$XX
    SSR.s= Y2.wait - pure.leaf[[i]]$Y2 - diag(t(xy.s) %*% ginv(xx.s) %*% xy.s)
    Dist.Gap = c(Dist.Gap, sum((SSR.all - SSR.s - pure.leaf[[i]]$SSR0)*weight))
    # Dist.Gap = c(Dist.Gap, sum((SSR.all - SSR.s - pure.leaf[[i]]$SSR0)/(SSR.s + pure.leaf[[i]]$SSR0)*(sum(obs.no) - 2*p.var)/p.var*weight))
    
    e.beta.s = e.beta - ginv(pure.leaf[[i]]$XX) %*% pure.leaf[[i]]$Xy 
    Dist.inter = c(Dist.inter, sum(t(e.beta.s) %*% pure.leaf[[i]]$XX %*% e.beta.s *weight))
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
    xy.s = Xy.wait - leaf.merge$Xy
    xx.s = XX.wait - leaf.merge$XX
    SSR.s= Y2.wait - leaf.merge$Y2 - diag(t(xy.s) %*% ginv(xx.s) %*% xy.s)
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
      for ( k in seq(nrow(Dist.tab)) ) {
        Y2.merge = pure.leaf[[k]]$Y2 + leaf.merge$Y2
        Xy.merge = pure.leaf[[k]]$Xy + leaf.merge$Xy
        XX.merge = pure.leaf[[k]]$XX + leaf.merge$XX
        # A.merge  = pure.leaf[[k]]$A.mat - pure.leaf[[k]]$A.mat %*% ginv(pure.leaf[[k]]$A.mat + leaf.merge$A.mat) %*% pure.leaf[[k]]$A.mat
        A.merge  = ginv(XX.merge)
        SSR.merge= diag(Y2.merge - t(Xy.merge) %*% A.merge %*% Xy.merge)
        Dist.new = c(Dist.new, sum(weight*((SSR.merge - pure.leaf[[k]]$SSR0 - leaf.merge$SSR0)/(pure.leaf[[k]]$SSR0 + leaf.merge$SSR0)*(pure.leaf[[k]]$N.in + leaf.merge$N.in - 2*p.var))))
      }
      Dist.tab = cbind(Dist.tab, Dist.new)
      Dist.tab = rbind(Dist.tab, rep(Inf, ncol(Dist.tab)))
    }
    
    if (parallel) {
      if (n.tmp <= stop) {
        break
      }
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
              p.var  = p.var, 
              e.sigma= sum(weight*e.sigma)))
}





