#' InstaPrism deconvolution with log-likelihood trace
# run devtools::load_all(path = 'path to InstaPrism package') to use this function

# For users interested in detailed model convergence status, we've introduced the `get_loglikelihood_trace()` function to visualize how log likelihood of the model converges over iterations. 
# `get_loglikelihood_trace()` function is the same as `InstaPrism()`` function except that it computes log likelihood of the model every  window_size iterations.
# example
# trace_obj = get_loglikelihood_trace(input_type = 'refPhi_cs',bulk_Expr = bulk_Expr,refPhi_cs = refPhi_obj,n.core = 16)
# plot(seq(1,length(trace_obj$loglikelihood_trace))*trace_obj$window_size,trace_obj$loglikelihood_trace)


Rcpp::sourceCpp('cppTraceFunc.cpp') # load bpFixedPointTrace() function

get_loglikelihood_trace = function(input_type=c('raw','prism','refPhi','refPhi_cs'),
                                   sc_Expr=NULL,bulk_Expr=NULL,
                                   cell.type.labels=NULL,cell.state.labels=NULL,
                                   filter=TRUE,
                                   outlier.cut=0.01,outlier.fraction=0.1,
                                   pseudo.min=1E-8,
                                   prismObj=NULL,
                                   refPhi=NULL,
                                   refPhi_cs=NULL,
                                   n.iter=NULL,
                                   window_size = 10,
                                   n.core=1){
  trace = function(bulk_Expr,ref,n.iter,n.core){
    cms=intersect(rownames(bulk_Expr),rownames(ref))
    ref=ref[cms,]
    pboptions(type = "txt", style = 3, char = "=")
    res.list=pbapply(bulk_Expr[cms,,drop=F],2,function(x) bpFixedPointTrace(bulk=matrix(x),ref=ref,n_iter = n.iter, window_size = window_size),cl = n.core)
    return(res.list)
  }
  
  
  if(input_type=='raw'){
    
    if(any(is.null(sc_Expr),is.null(cell.type.labels),is.null(cell.state.labels),is.null(bulk_Expr))){
      stop('Need to specify all raw input objects. One or more input objects from the following is missing: sc_Expr, cell.type.labels, cell.state.labels, bulk_Expr')
    }
    
    if(length(commonRows(sc_Expr,bulk_Expr)) < 10){
      stop('few gene overlap detected between sc_Expr and bulk_Expr, please ensure consistent gene symbol formats')
    }
    
    if(is.null(n.iter)){
      n.iter = max(100, 2*length(unique(cell.state.labels)))
    }
    
    
    bp=bpPrepare(input_type='raw',
                 sc_Expr,
                 cell.type.labels,cell.state.labels,
                 bulk_Expr,filter,
                 outlier.cut,outlier.fraction,
                 pseudo.min)
    
    cat('deconvolution with scRNA reference phi \n')
    
    res.list = trace(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)
    cell.states = colnames(bp@phi.cs)
    map=bp@map
  }else if(input_type=='prism'){
    
    if(is.null(prismObj)){
      stop('Need to specify a prismObj when input_type = "prism"')
    }
    
    if(is.null(n.iter)){
      n.iter = max(100, 2*nrow(prismObj@phi_cellState@phi))
    }
    
    cat('deconvolution with scRNA reference phi \n')
    res.list = trace(t(prismObj@mixture),t(prismObj@phi_cellState@phi),n.iter,n.core)
    cell.states = rownames(prismObj@phi_cellState@phi)
    map=prismObj@map
  }else if(input_type=='refPhi'){
    if(any(is.null(refPhi),is.null(bulk_Expr))){
      stop('Need to specify refPhi and bulk_Expr when input_type = "refPhi"')
    }
    
    if(length(commonRows(refPhi@phi.ct,bulk_Expr)) < 10){
      stop('few gene overlap detected between reference and bulk_Expr, please ensure consistent gene symbol formats')
    }
    
    if(is.null(n.iter)){
      n.iter = max(100, 2* ncol(refPhi@phi.cs))
    }
    
    bp = bpPrepare(input_type = 'refPhi',
                   bulk_Expr = bulk_Expr,
                   filter = filter,
                   outlier.cut=outlier.cut,
                   outlier.fraction=outlier.fraction,
                   pseudo.min=pseudo.min,
                   refPhi = refPhi)
    
    cat('deconvolution with scRNA reference phi \n')
    res.list = trace(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)
    cell.states = colnames(bp@phi.cs)
    map=bp@map
    
  }else if(input_type=='refPhi_cs'){
    if(any(is.null(refPhi_cs),is.null(bulk_Expr))){
      stop('Need to specify refPhi_cs and bulk_Expr when input_type = "refPhi_cs"')
    }
    
    if(length(commonRows(refPhi_cs@phi.cs,bulk_Expr)) < 10){
      stop('few gene overlap detected between reference and bulk_Expr, please ensure consistent gene symbol formats')
    }
    
    if(is.null(n.iter)){
      n.iter = max(100, 2* ncol(refPhi_cs@phi.cs))
    }
    
    
    bp = bpPrepare(input_type = 'refPhi_cs',
                   bulk_Expr = bulk_Expr,
                   filter = filter,
                   outlier.cut=outlier.cut,
                   outlier.fraction=outlier.fraction,
                   pseudo.min=pseudo.min,
                   refPhi_cs = refPhi_cs)
    
    cat('deconvolution with scRNA reference phi \n')
    res.list = trace(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)
    cell.states = colnames(bp@phi.cs)
    map=bp@map
  }
  
  theta=mapply(`[[`, res.list, 1)
  rownames(theta)=cell.states
  theta.ct = do.call(rbind,lapply(map,function(x)colSums(theta[rownames(theta) %in% x,,drop=F])))
  
  loglikelihood_trace = rowSums(mapply(`[[`, res.list, 2))
  
  return(list(Post.ini.cs = theta,
              Post.ini.ct = theta.ct,
              loglikelihood_trace = loglikelihood_trace,
              n.iter = n.iter,
              window_size = window_size))
  
}