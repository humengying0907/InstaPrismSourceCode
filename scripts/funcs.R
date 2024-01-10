set_rowNames_to_NULL = function(x){
  rownames(x) = NULL
  return(x)
}


UMI2CPM = function(scExpr){
  col_sums <- Matrix::colSums(scExpr)
  scaling_factors <- 10^6 / col_sums
  cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
  colnames(cpm) = colnames(scExpr)
  return(cpm)
}

factor2vector_dataframe = function(meta){
  factor_columns <- sapply(meta, class) == "factor"
  meta[factor_columns] <- lapply(meta[factor_columns], as.character)
  return(meta)
}


unlog_transform <- function(x,base = 2) {
  if (x == 0) {
    return(0)  # If the element is 0, return 0
  } else {
    return(base^x-1)  # Apply the inverse transformation
  }
}

# using refPhi
multi_deconv = function(bulk_simulation_obj, refPhi_cs_obj, key, updateReference,niter,n.core){
  
  deconvResults = list()
  for(bulk in names(bulk_simulation_obj)){
    bulk_expr = bulk_simulation_obj[[bulk]]$simulated_bulk
    bulk_res = list()
    bulk_res[['initial']] = InstaPrism(input_type = 'refPhi_cs',
                                       bulk_Expr = bulk_expr,
                                       refPhi_cs = refPhi_cs_obj,
                                       n.iter = niter,
                                       n.core = 16)
    if(updateReference){
      bulk_res[['updated']] = InstaPrism_update(bulk_res[['initial']],
                                                bulk_Expr = bulk_expr,
                                                n.iter = niter,
                                                n.core = n.core,
                                                cell.types.to.update = 'all',
                                                key = key)
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
    bulk_theta_list[['initial']] = InstaPrism_res_obj[[bulk]]$initial@Post.ini.ct@theta %>% t()
    
    updateModes = names(InstaPrism_res_obj[[bulk]])[-1]
    for(mode in updateModes){
      updated_obj = InstaPrism_res_obj[[bulk]][[mode]]
      
      # determine if updated_obj@theta contain cell.state information
      
      if(all(rownames(updated_obj@theta) %in% colnames(bulk_theta_list[['initial']]))){
        bulk_theta_list[[mode]] = updated_obj@theta %>% t()
      }else{
        bulk_theta_list[[mode]] = merge_updated_theta(updated_obj) %>% t()
      }
    }
    theta_list[[bulk]] = bulk_theta_list
  }
  return(theta_list)
}

theta_evalu = function(theta_list,sim_obj,nonzero_threshold = 0.05){
  deconvPerformance = list()
  
  for(bulk in names(theta_list)){
    bulk_performance = list()
    
    Y = sim_obj[[bulk]]$simulated_frac
    Y = Y[,order(colnames(Y)),drop = F]
    
    # Only cell types with a non-zero count ratio equal to or higher than 'nonzero_ratios' will be included for comparing their deconvolution performance.
    nonzero_ratios = colSums(Y>0, na.rm = T)/nrow(Y)
    ct_to_evalu = names(nonzero_ratios)[nonzero_ratios >= nonzero_threshold]
    Y_ = Y[,ct_to_evalu,drop = F]

    # detailed per-method evaluation statistics
    deconvEvalu = list()
    
    # look at the distribution of predicted theta for cell types that are in the reference, but not in bulk simulation
    over_representation = list()
    
    for(method in names(theta_list[[bulk]])){
      E = theta_list[[bulk]][[method]]
      E = E[,order(colnames(E)),drop = F]
      
      if(ncol(Y) == ncol(E)){
        if(all(colnames(Y) == colnames(E))){
          # modify E according to nonzero_threshold if it contains the exact set of cell types as the true_frac
          E = E[,ct_to_evalu,drop = F]
          
          op_df = data.frame(matrix(NA,ncol = 2,nrow = 0))
          colnames(op_df) = c('over_represented_ct','estimate')
          over_representation[[method]] = op_df
          
        }else{
          # find cell types that exhibit the highest correlation with true_frac in Y_
          maxCorName = c()
          for(ct in colnames(Y_)){
            id = apply(cor(E,Y_[,ct],use = 'pairwise.complete.obs'),2,which.max)
            maxCorName = c(maxCorName, colnames(E)[id])
          }
          
          op_ct = colnames(E)[!colnames(E) %in% unique(maxCorName)]
          op_df = gather(E[,op_ct,drop=F] %>% as.data.frame(),over_represented_ct,estimate)
          over_representation[[method]] = op_df
          
          E = E[,maxCorName,drop = F]
          colnames(E) = make.unique(colnames(E))
          
        }
      }else{
        maxCorName = c()
        for(ct in colnames(Y_)){
          id = apply(cor(E,Y_[,ct],use = 'pairwise.complete.obs'),2,which.max)
          maxCorName = c(maxCorName, colnames(E)[id])
        }
        


        op_ct = colnames(E)[!colnames(E) %in% unique(maxCorName)]
        op_df = gather(E[,op_ct,drop=F] %>% as.data.frame(),over_represented_ct,estimate)
        over_representation[[method]] = op_df
        
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
    
    
    bulk_performance$deconvEvalu = deconvEvalu
    bulk_performance$over_representation = over_representation
    
    deconvPerformance[[bulk]] = bulk_performance
    message(paste('>>>>>>>>>>>>>>>>>>>>>>>>>>>>> finish performance evaluation for',bulk,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))
    
  }
  return(deconvPerformance)
  
}

theta_plot = function(single_theta_performance,method){
  l = single_theta_performance$deconvEvalu[[method]]
  op_df = single_theta_performance$over_representation[[method]]
  op_df$kind = 'ct in reference, not in bulk'
  
  M=l$M
  summ=l$summ
  summ = summ %>% rownames_to_column('cell_type')
  nrow = 1
  
  title = method
  
  p = ggplot(M,aes(true_frac,estimate))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1,linetype='dotted',color='red')+
    facet_wrap(~cell_type,nrow=nrow)+  # scatter plot will be arranged alphabetically
    theme(aspect.ratio=1)+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  p + ggpp::geom_table_npc(data = summ, label = lapply(split(summ, summ$cell_type),
                                                            FUN = function(entry) {subset(entry, select = -cell_type)}),
                                npcx = 0.00, npcy = 1, hjust = 0, vjust = 1, size=3,
                                table.theme = ttheme_gtlight)
}

theta_op_plot = function(single_theta_performance,method){
  l = single_theta_performance$deconvEvalu[[method]]
  op_df = single_theta_performance$over_representation[[method]]
  op_df$kind = 'ct in reference, not in bulk'
  
  M=l$M
  summ=l$summ
  summ = summ %>% rownames_to_column('cell_type')
  nrow = 1
  
  ggplot(op_df,aes(over_represented_ct,estimate))+
    geom_boxplot()+
    ylim(c(0,1))+
    facet_wrap(~kind)+
    theme_bw()
  }

theta_plot_export = function(single_theta_performance,fig_width = 12,fig_height = 18,plot_file_name = 'theta_performance.pdf'){
  
  plot_sets = list()
  
  for(method in names(single_theta_performance$deconvEvalu)){
    
    l = single_theta_performance$deconvEvalu[[method]]
    #op_df = single_theta_performance$over_representation[[method]]
    #op_df$kind = 'ct in reference, not in bulk'
    
    M=l$M
    summ=l$summ
    summ = summ %>% rownames_to_column('cell_type')
    nrow = 1
    
    title = method
    
    #if(grepl('updateMode',method)){
    #  mode_id = as.numeric(gsub("[^0-9]", "", method))
    #  title = paste0(method,' : ',key,', ',paste(ct.update.list[[mode_id]],collapse = ', '))
    #}else{
    #  title = method
    #}
    
    p = ggplot(M,aes(true_frac,estimate))+
      geom_point()+
      geom_abline(intercept = 0, slope = 1,linetype='dotted',color='red')+
      facet_wrap(~cell_type,nrow=nrow)+  # scatter plot will be arranged alphabetically
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


create_pseudobulk_obj = function(scExpr,scMeta,colnames_of_sample = 'sample',colnames_of_cellType = 'cell_type',
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




