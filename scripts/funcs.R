factor2vector_dataframe = function(meta){
  factor_columns <- sapply(meta, class) == "factor"
  meta[factor_columns] <- lapply(meta[factor_columns], as.character)
  return(meta)
}

map_count = function(cell.type.labels,cell.state.labels){
  df = data.frame(cell_type = cell.type.labels,
                  cell_state = cell.state.labels)
  v = df %>% group_by(cell_type,cell_state) %>% summarise(n=n())
  return(as.data.frame(v))
}

# using refPhi
multi_deconv = function(bulk_simulation_obj, refPhi_list, key, updateReference,niter,n.core){
  
  deconvResults = list()
  for(bulk in names(bulk_simulation_obj)){
    bulk_expr = bulk_simulation_obj[[bulk]]$simulated_bulk
    bulk_res = list()
    
    for(refName in names(refPhi_list)){
      refPhi_cs_obj = refPhi_list[[refName]]
      bulk_res[[paste0(refName)]] = InstaPrism(input_type = 'refPhi_cs',
                                                           bulk_Expr = bulk_expr,
                                                           refPhi_cs = refPhi_cs_obj,
                                                           n.iter = niter,
                                                           n.core = n.core)
      if(updateReference){
        bulk_res[[paste0(refName,'_','updated')]] = InstaPrism_update(bulk_res[[paste0(refName)]],
                                                  bulk_Expr = bulk_expr,
                                                  n.iter = niter,
                                                  n.core = n.core,
                                                  cell.types.to.update = 'all',
                                                  key = key)
      }
    }
    
    deconvResults[[bulk]] = bulk_res
    
    print(paste('+++++++++++++++++++++++++++++++++++++++++++++ finish',bulk,'+++++++++++++++++++++++++++++++++++++++++++++'))
  }
  return(deconvResults)
  
}

extract_theta = function(InstaPrism_res_obj){
  theta_list = list()
  for (bulk in names(InstaPrism_res_obj)){
    bulk_theta_list = list()
    for(refName in names(InstaPrism_res_obj[[bulk]])){
      if(grepl('_updated',refName)){
        bulk_theta_list[[refName]] = InstaPrism_res_obj[[bulk]][[refName]]@theta %>% t()
      }else{
        bulk_theta_list[[refName]] = InstaPrism_res_obj[[bulk]][[refName]]@Post.ini.ct@theta %>% t()
      }
    }
    
    theta_list[[bulk]] = bulk_theta_list
  }
  return(theta_list)
}

theta_evalu = function(theta_list,simulated_frac,nonzero_threshold = 0.05){
  
  Y = simulated_frac
  Y = Y[,order(colnames(Y)),drop = F]
  
  # Only cell types with a non-zero count ratio equal to or higher than 'nonzero_ratios' will be included for comparing their deconvolution performance.
  nonzero_ratios = colSums(Y>0.005, na.rm = T)/nrow(Y)
  ct_to_evalu = names(nonzero_ratios)[nonzero_ratios >= nonzero_threshold]
  Y_ = Y[,ct_to_evalu,drop = F]
  
  # detailed per-method evaluation statistics
  deconvEvalu = list()
  
  for(method in names(theta_list)){
    
    E = theta_list[[method]]
    E = E[,order(colnames(E)),drop = F]
    
    if(ncol(Y) == ncol(E)){
      if(all(colnames(Y) == colnames(E))){
        # modify E according to nonzero_threshold if it contains the exact set of cell types as the true_frac
        E = E[,ct_to_evalu,drop = F]
        
      }else{
        # find cell types that exhibit the highest correlation with true_frac in Y_
        maxCorName = c()
        for(ct in colnames(Y_)){
          id = apply(cor(E,Y_[,ct],use = 'pairwise.complete.obs'),2,which.max)
          maxCorName = c(maxCorName, colnames(E)[id])
        }
        
        E = E[,maxCorName,drop = F]
        colnames(E) = make.unique(colnames(E))
        
      }
    }else{
      maxCorName = c()
      for(ct in colnames(Y_)){
        id = apply(cor(E,Y_[,ct],use = 'pairwise.complete.obs'),2,which.max)
        maxCorName = c(maxCorName, colnames(E)[id])
      }
      
      E = E[,maxCorName,drop = F]
      colnames(E) = make.unique(colnames(E))
    }
    
    m1 = gather(Y_ %>% as.data.frame(),cell_type,true_frac)
    m2 = gather(E %>% as.data.frame(),maxCorName,estimate)
    M = cbind(m1,m2)
    
    # add summary statistics
    summ <- M %>%
      group_by(cell_type) %>%
      summarise(
        RMSE = caret::RMSE(true_frac, estimate,na.rm = T),
        cor = cor(true_frac,estimate,method = 'pearson',use = 'pairwise.complete.obs')) %>%
      mutate_if(is.numeric, round, digits=2) %>% as.data.frame() %>% column_to_rownames('cell_type')
    
    summ$maxCorName = M$maxCorName[match(rownames(summ),M$cell_type)]
    deconvEvalu[[method]] = list(M = M, summ = summ)
  }
  
  return(deconvEvalu)
  
}

theta_plot_export = function(single_theta_performance,fig_width = 12,fig_height = 18,
                             plot_file_name = 'theta_performance.pdf'){
  
  plot_sets = list()
  
  for(method in names(single_theta_performance)){
    
    l = single_theta_performance[[method]]
    M=l$M
    summ=l$summ
    summ = summ %>% rownames_to_column('cell_type')
    nrow = 1
    summ = summ[order(summ$cell_type),] # order alphabetically

    title = method
    
    p = ggplot(M,aes(true_frac,estimate))+
      geom_point()+
      geom_abline(intercept = 0, slope = 1,linetype='dotted',color='red')+
      facet_wrap(~cell_type,nrow=nrow)+  
      ggtitle(title)+
      ylim(c(0,1.6))+
      theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1)+
      theme_bw()
    
    p0 = p + ggpp::geom_table_npc(data = summ, label = lapply(split(summ, summ$cell_type),
                                                              FUN = function(entry) {subset(entry, select = -cell_type)}),
                                  npcx = 0.00, npcy = 1, hjust = 0, vjust = 1, size=3,
                                  table.theme = ttheme_gtlight)
    
    plot_sets = c(plot_sets,list(p0))
    
  }
  p_out=gridExtra::grid.arrange(grobs = lapply(
    plot_sets,
    egg::set_panel_size,
    width = unit(4.5, "cm"),
    height = unit(3.75, "cm")),ncol = 1
  )
  
  ggsave(plot_file_name,width = fig_width,height = fig_height,plot = p_out)
  
}

create_pseudobulk_obj = function(scExpr,scMeta,colnames_of_sample = 'sample',
                                 colnames_of_cellType = 'cell_type',
                                 unit = 'cpm',min.cells = 3,
                                 add_heter = F, simulated_frac = NULL,min_chunkSize = 10){
  
  # keep genes that have expression values in at least min.cells
  v = Matrix::rowSums(scExpr > 0)
  scExpr = scExpr[v >= min.cells,]
  
  # keep samples with cell type assignment
  to_keep = which(scMeta[,colnames_of_cellType] != '')
  scExpr = scExpr[,to_keep]
  scMeta = scMeta[to_keep,]
  
  cell_usage = scMeta[,colnames_of_sample,drop = F] %>% rownames_to_column('cell')
  colnames(cell_usage)[2] = 'bulk_id'
  sample_identifiers = scMeta[,colnames_of_sample]
  sample_identifiers = as.character(sample_identifiers)
  
  # use cell count as approximation to fractions
  Meta = data.frame(row.names = rownames(scMeta),
                    cell_type = scMeta[,colnames_of_cellType],
                    sampleID = scMeta[,colnames_of_sample])
  
  cell_usage$cell_type = Meta$cell_type[match(cell_usage$cell,rownames(Meta))]
  
  
  cell_type_table = Meta %>%
    group_by(sampleID, cell_type) %>%
    summarise(n=n()) %>%
    mutate(frac=n/sum(n)) %>% as.data.frame()
  frac = spread(cell_type_table[,-3], cell_type, frac) %>% column_to_rownames('sampleID') %>% as.matrix()
  frac[is.na(frac)] = 0
  
  if(unit == 'cpm'){
    # rowMeans of cells from the same sample
    pseudobulk = deconvBenchmarking:::build_ref_matrix(scExpr,sample_identifiers)
    pseudobulk = pseudobulk[,rownames(frac)]
    
    if(add_heter){
      heter = bulkSimulator_heter(scExpr = scExpr,
                                  scMeta = Meta,
                                  colnames_of_cellType = 'cell_type',
                                  colnames_of_sample = 'sampleID',
                                  simulated_frac = simulated_frac,
                                  use_chunk = 'random',
                                  export_cellUsage = T,
                                  n.core = 16)
      pseudobulk_obj = list()
      pseudobulk_obj$simulated_bulk = cbind(heter$simulated_bulk,pseudobulk)
      pseudobulk_obj$simulated_frac = rbind(heter$simulated_frac,frac)
      pseudobulk_obj$cell_usage = rbind(heter$cell_usage,cell_usage)
    }else{
      pseudobulk_obj = list(simulated_bulk = pseudobulk,
                            simulated_frac = frac,
                            cell_usage = cell_usage)
    }
    
  }else if(unit == 'UMI'){
    col_sums <- Matrix::colSums(scExpr)
    scaling_factors <- 10^6 / col_sums
    cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
    colnames(cpm) = colnames(scExpr)
    pseudobulk = deconvBenchmarking:::build_ref_matrix(cpm,sample_identifiers)
    pseudobulk = pseudobulk[,rownames(frac)]
  
    
    if(add_heter){
      heter = bulkSimulator_heter(scExpr = cpm,
                                  scMeta = Meta,
                                  colnames_of_cellType = 'cell_type',
                                  colnames_of_sample = 'sampleID',
                                  simulated_frac = simulated_frac,
                                  min_chunkSize = min_chunkSize,
                                  use_chunk = 'random',
                                  export_cellUsage = T,
                                  n.core = 16)
      pseudobulk_obj = list()
      pseudobulk_obj$simulated_bulk = cbind(heter$simulated_bulk,pseudobulk)
      pseudobulk_obj$simulated_frac = rbind(heter$simulated_frac,frac)
      pseudobulk_obj$cell_usage = rbind(heter$cell_usage,cell_usage)
      
    }else{
      pseudobulk_obj= list(simulated_bulk = pseudobulk,
                           simulated_frac = frac,
                           cell_usage = cell_usage)
    }
    
  }
  
  # raise warning if medians of simulated_bulk is close to 0: indicating not enough cells to create the pseudobulk objects
  col_medians = apply(pseudobulk_obj$simulated_bulk,2,median)
  if(any(col_medians<0.1)){
    poor_quality_count = sum(col_medians < 0.1)
    warning(paste(poor_quality_count,'samples have expression median less than 0.1'))
  }
  
  return(pseudobulk_obj)
}




