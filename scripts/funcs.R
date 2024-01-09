set_rowNames_to_NULL = function(x){
  rownames(x) = NULL
  return(x)
}

build_ref_matrix2 = function(Expr,cell_type_labels){
  stopifnot(ncol(Expr)==length(cell_type_labels))
  group = list()
  for(i in unique(cell_type_labels)){
    group[[i]] <- which(cell_type_labels %in% i)
  }
  C = do.call(cbind, lapply(group,function(x) Matrix::rowSums(Expr[,x,drop=F])))
  C
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


quick_map = function(cell.type.labels,cell.state.labels){
  cell.type.labels = as.vector(cell.type.labels)
  cell.state.labels = as.vector(cell.state.labels)
  map=list()
  for (ct in unique(cell.type.labels)){
    i=which(cell.type.labels==ct)
    map[[ct]]=unique(cell.state.labels[i])
  }
  return(map)
}

map_count = function(cell.type.labels,cell.state.labels){
  df = data.frame(cell_type = cell.type.labels,
                  cell_state = cell.state.labels)
  v = df %>% group_by(cell_type,cell_state) %>% summarise(n=n())
  return(as.data.frame(v))
}


map_alluvial = function(cell.type.labels,cell.state.labels){
  require(ggplot2)
  require(ggalluvial)
  df = data.frame(cell_type = cell.type.labels,
                  cell_state = cell.state.labels)
  
  data_long <- as.data.frame(table(df))
  ggplot(data_long,
         aes(axis1 = cell_type, axis2 = cell_state, y = Freq)) +
    geom_alluvium(aes(fill = cell_type)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    theme(legend.position = 'none')
}

meta_summary = function(scMeta,colnames_of_sample,colnames_of_cellType){
  scMeta_summary=scMeta %>%   group_by_(colnames_of_sample, colnames_of_cellType) %>%
    summarise(n=n()) %>%
    mutate(freq=n/sum(n))
  return(scMeta_summary)
}

refPhi_heatmap = function(phi){
  require(ggcorrplot)
  norm_factor = rowSums(phi)
  pseudo.min = min(phi[phi>0])
  norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)
  phi_rownormalized = sweep(phi,1,norm_factor,'/')
  ggcorrplot(cor(phi_rownormalized),lab = T,lab_size = 3)
}

ggboxplot = function(df){
  df = df %>% as.data.frame()
  df_long <- gather(df, key = "Variable", value = "Value")
  ggplot(df_long, aes(x = reorder(Variable,Value,median), y = Value)) +
    geom_boxplot()+
    coord_flip()+
    theme_bw()
}


unlog_transform <- function(x,base = 2) {
  if (x == 0) {
    return(0)  # If the element is 0, return 0
  } else {
    return(base^x-1)  # Apply the inverse transformation
  }
}

# using refPhi
get_InstaPrism_res = function(bulk_simulation_obj, refPhi, key, ct.update.list){

  deconvResults = list()
  for(bulk in names(bulk_simulation_obj)){
    bulk_expr = bulk_simulation_obj[[bulk]]$simulated_bulk
    bulk_res = list()
    bulk_res[['initial']] = InstaPrism(input_type = 'refPhi',
                                       bulk_Expr = bulk_expr,
                                       refPhi = refPhi,
                                       n.iter = 300,
                                       n.core = 16)
    
    for(i in 1:length(ct.update.list)){
      bulk_res[[paste0('updateMode_',i)]]=InstaPrism_update(bulk_res[['initial']],
                                                            bulk_Expr = bulk_expr,
                                                            n.iter = 300,
                                                            n.core = 16,
                                                            cell.types.to.update = ct.update.list[[i]],
                                                            key = key, 
                                                            keep.phi = 'phi.cs')
    }
    
    deconvResults[[bulk]] = bulk_res
    
    print(paste('+++++++++++++++++++++++++++++++++++++++++++++ finish',bulk,'+++++++++++++++++++++++++++++++++++++++++++++'))
  }
  return(deconvResults)
  
}

# using refPhi_cs (for refPhi_cs validataion)
multi_deconv = function(bulk_simulation_obj, refPhi_cs_obj, key, ct.update.list,niter = 300){
  
  deconvResults = list()
  for(bulk in names(bulk_simulation_obj)){
    bulk_expr = bulk_simulation_obj[[bulk]]$simulated_bulk
    bulk_res = list()
    bulk_res[['initial']] = InstaPrism(input_type = 'refPhi_cs',
                                       bulk_Expr = bulk_expr,
                                       refPhi_cs = refPhi_cs_obj,
                                       n.iter = niter,
                                       n.core = 16)
    
    for(i in 1:length(ct.update.list)){
      bulk_res[[paste0('updateMode_',i)]]=InstaPrism_update(bulk_res[['initial']],
                                                            bulk_Expr = bulk_expr,
                                                            n.iter = niter,
                                                            n.core = 16,
                                                            cell.types.to.update = ct.update.list[[i]],
                                                            key = key, 
                                                            keep.phi = 'phi.cs')
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



reference_comparison = function(InstaPrism_res_obj,
                                theta_performance_obj,
                                ground_truth_refPhi, 
                                bulk_ct_env,
                                ct.update.list){
  ref_comparison_list = list()
  for(bulk in names(InstaPrism_res_obj)){
    methods = names(InstaPrism_res_obj[[bulk]])
    
    # first consider bk_ct (cell types that actually present in bulk)
    bk_ct_phi_comparison = data.frame(matrix(NA,ncol = 7, nrow = 0))
    colnames(bk_ct_phi_comparison) = c('gene','phi_true','phi_true_adjusted','phi_inferred','cell_type','maxCorName','method')
    
    for(ct in bulk_ct_env){
      print(ct)
      
      ct_phi_true = ground_truth_refPhi@phi.ct[,ct]
      
      get_ct_ref = function(method){
        
        summ = theta_performance_obj[[bulk]]$deconvEvalu[[method]]$summ
        maxCorName = summ$maxCorName[which(rownames(summ)==ct)]
        
        
        if(grepl('\\d',method)){
          mode_id = as.numeric(gsub("[^0-9]", "", method))
          updated_ct = ct.update.list[[mode_id]]
          if(maxCorName %in% updated_ct){
            ct_phi_ref =  InstaPrism_res_obj[[bulk]][[method]]@psi_env[,maxCorName]
            cms = intersect(names(ct_phi_true),names(ct_phi_ref))
            
            ct_phi_true = ct_phi_true[cms]
            ct_phi_true_adjusted = ct_phi_true/sum(ct_phi_true)
            ct_phi_ref = ct_phi_ref[cms]
            
            
            out = data.frame(phi_true = ct_phi_true,
                             phi_true_adjusted = ct_phi_true_adjusted,
                             phi_inferred = ct_phi_ref) %>% rownames_to_column('gene') %>% mutate(cell_type = ct, maxCorName = maxCorName, method = method)
            
          }else{
            out = data.frame(matrix(NA,ncol = 7, nrow = 0))
            colnames(out) = c('gene','phi_true','phi_true_adjusted','phi_inferred','cell_type','maxCorName','method')
          }
        }else{
          ct_phi_ref = InstaPrism_res_obj[[bulk]][[method]]@initial.reference@phi.ct[,maxCorName]
          
          
          cms = intersect(names(ct_phi_true),names(ct_phi_ref))
          
          ct_phi_true = ct_phi_true[cms]
          ct_phi_true_adjusted = ct_phi_true/sum(ct_phi_true)
          ct_phi_ref = ct_phi_ref[cms]
          
          out = data.frame(phi_true = ct_phi_true,
                           phi_true_adjusted = ct_phi_true_adjusted,
                           phi_inferred = ct_phi_ref) %>% rownames_to_column('gene') %>% mutate(cell_type = ct, maxCorName = maxCorName, method = method)
        }
        return(out)
      }
      
      ct_phi_comparison = do.call(rbind,lapply(methods,get_ct_ref))
      bk_ct_phi_comparison = rbind(bk_ct_phi_comparison,ct_phi_comparison)
      
    }
    
    bk_ct_phi_comparison$log_phi_true = log2(bk_ct_phi_comparison$phi_true + 1e-08)
    bk_ct_phi_comparison$log_phi_true_adjusted = log2(bk_ct_phi_comparison$phi_true_adjusted + 1e-08)
    bk_ct_phi_comparison$log_phi_inferred = log2(bk_ct_phi_comparison$phi_inferred + 1e-08)
    
    
    ref_comparison_list[[bulk]]$bk_ct_phi_comparison = bk_ct_phi_comparison
    
    # next consider ct over-represented in refPhi: consider the simple scenario that op_ct is the same across update strategies
    ref_ct = names(InstaPrism_res_obj[[bulk]]$initial@map)
    summ_initial = theta_performance_obj[[bulk]]$deconvEvalu$initial$summ
    op_ct = ref_ct[!ref_ct %in% summ_initial$maxCorName]
    
    op_ct_initial = InstaPrism_res_obj[[bulk]]$initial@initial.reference@phi.ct[,op_ct]
    genes = rownames(op_ct_initial)
    op_ct_initial = op_ct_initial %>% as.data.frame() %>% rownames_to_column('gene') %>% gather(cell_type,initial_phi,2:(length(op_ct)+1))
    
    methods_update = methods[grepl('\\d',methods)]
    op_ct_phi_updated_list = list()
    for(method in methods_update){
      print(method)
      
      op_ct_phi = InstaPrism_res_obj[[bulk]][[method]]@psi_env[,op_ct]
      op_ct_phi = op_ct_phi[genes,]
      op_ct_phi = op_ct_phi %>% as.data.frame() %>% rownames_to_column('gene') %>% gather(cell_type,updated_phi,2:(length(op_ct)+1))
      stopifnot(all.equal(op_ct_phi$gene,op_ct_initial$gene))
      stopifnot(all.equal(op_ct_phi$cell_type,op_ct_initial$cell_type))
      
      op_ct_phi$initial_phi = op_ct_initial$initial_phi
      op_ct_phi$method = method
      
      op_ct_phi_updated_list[[method]] = op_ct_phi
    }
    
    op_ct_phi_comparison = do.call(rbind,op_ct_phi_updated_list) 
    rownames(op_ct_phi_comparison) = NULL
    op_ct_phi_comparison$log_initial_phi = log2(op_ct_phi_comparison$initial_phi+1e-08)
    op_ct_phi_comparison$log_updated_phi = log2(op_ct_phi_comparison$updated_phi+1e-08)
    
    
    ref_comparison_list[[bulk]]$op_ct_phi_comparison = op_ct_phi_comparison
    
    message(paste('>>>>>>>>>>>>>>>>>>>>>>>>>>>>> finish reference comparison for',bulk,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))
    
  }
  return(ref_comparison_list)
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

Z_ground_truth = function(scExpr, 
                          simulated_bulk,
                          cell_usage,
                          unit = 'cpm',
                          n.core = 16){
  
  scExpr = scExpr[rownames(simulated_bulk),]
  bulk_ids = colnames(simulated_bulk)
  
  Z_truth = list()
  
  if(unit == 'UMI'){
    col_sums <- Matrix::colSums(scExpr)
    scaling_factors <- 10^6 / col_sums
    cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
    colnames(cpm) = colnames(scExpr)
    
    scExpr = cpm
  }
  
  Z_ct_single = function(bulk_id,scExpr,cell_usage_sub){
    cells = cell_usage_sub$cell[cell_usage_sub$bulk_id == bulk_id]
    if(length(cells)==0){
      avg_expr = rep(0,nrow(scExpr))
      names(avg_expr) = rownames(scExpr)
      return(avg_expr)
    }else{
      avg_expr = rowMeans(scExpr[,cells,drop=F])
      return(avg_expr)
    }
  }
  
  Z_ct_multi = function(bulk_ids,scExpr,cell_usage_sub,n.core){
    
    mat = do.call(cbind,pblapply(bulk_ids,Z_ct_single,scExpr,cell_usage_sub,cl = n.core))
    colnames(mat) = bulk_ids
    return(mat)
  }
  
  bulk_ct_components = unique(cell_usage$cell_type)
  for(ct in bulk_ct_components){
    print(ct)
    cell_usage_sub = cell_usage[cell_usage$cell_type == ct,]
    
    Z_truth[[ct]] = Z_ct_multi(bulk_ids,scExpr,cell_usage_sub,n.core)
  }
  
  return(Z_truth)
  
}


pair_cor = function(a,b,margin = 'row',cor.method = 'pearson'){
  # a: gene * sample
  # b: gene * sample
  
  # make sure a,b has the same rows and columns
  stopifnot(all.equal(colnames(a),colnames(b)))
  stopifnot(all.equal(rownames(a),rownames(b)))
  
  
  corr = c()
  
  if(margin == 'row'){
    for(i in 1:nrow(a)){
      corr[i] = cor(a[i,],b[i,],method = cor.method,use = 'pairwise.complete.obs')
    }
  }else if(margin == 'column'){
    for(i in 1:ncol(a)){
      corr[i] = cor(a[,i],b[,i],method = cor.method,use = 'pairwise.complete.obs')
      
    }
  }

  
  return(corr)
}

read.gmt<-function(gmt.file){
  g<-GSA::GSA.read.gmt(gmt.file)
  genelist<-g$genesets
  names(genelist)<-g$geneset.names
  rm(g)
  return(genelist)
}

singscore_result=function(Expr,gsva_list){
  # filter element in the list with at least 2 elements
  filter_criteria= do.call(c,lapply(gsva_list,length))
  gsva_list=gsva_list[names(filter_criteria)[filter_criteria>=2]]
  
  rankData<-rankGenes(Expr)
  quick_singscore=function(genes,rankData){
    df<-simpleScore(rankData, upSet =genes)
    df$TotalScore
  }
  m=do.call(rbind,lapply(gsva_list,quick_singscore,rankData))
  colnames(m)=colnames(Expr)
  return(m)
}



Z_evalu = function(InstaPrism_res_obj,
                   theta_performance_obj,
                   sim_obj,
                   Z_truth_obj,
                   halmark.genelist,
                   MP_programs,
                   MP_map,
                   n.core = 16){
  bulk_names = names(InstaPrism_res_obj)
  bulk_names = bulk_names[bulk_names!='TCGA']
  
  Z_evalu_obj = list()
  
  hallmark_genes = do.call(c,halmark.genelist) %>% unique() # only compare within-gene variance for these genes
  
  
  for(bulk in bulk_names){
    simulated_frac = sim_obj[[bulk]]$simulated_frac
    simulated_bulk = sim_obj[[bulk]]$simulated_bulk
    
    summ = theta_performance_obj[[bulk]]$deconvEvalu$initial$summ # only consider initial for now
    summ$maxCorName = sub("\\.[^.]*$", "", summ$maxCorName)
    
    # if multiple ct match to the same reference ct, keep the one with the highest cor
    if(sum(duplicated(summ$maxCorName))>0){
      duplicate_indices <- which(duplicated(summ$maxCorName) | duplicated(summ$maxCorName, fromLast = TRUE))
      
      summ_a = summ[-duplicate_indices,]
      summ_b = summ[duplicate_indices,]
      
      summ_b = summ_b %>% rownames_to_column('ground_truth_ct') %>%
        group_by(maxCorName) %>%
        slice(which.max(cor)) %>%
        ungroup()  %>% column_to_rownames('ground_truth_ct')
      
      summ = rbind(summ_a,summ_b )
    }
    
    between_sample_cor = data.frame(matrix(NA,nrow = 0,ncol = 4))
    between_sample_MP_cor = data.frame(matrix(NA,nrow = 0,ncol = 4))
    
    within_sample_cor = data.frame(matrix(NA,nrow = 0,ncol = 4))
    
    
    for(ct in rownames(summ)){
      print(ct)
      maxCorName = summ$maxCorName[rownames(summ) == ct]
      Z_truth = Z_truth_obj[[bulk]][[ct]]
      Z_inferred = reconstruct_Z_ct_initial(InstaPrism_res_obj[[bulk]]$initial,maxCorName,n.core)
      
      # only use samples with sim-frac > 0 (otherwise Z-truth will be all zeros)
      to_use_samples = which(simulated_frac[,ct]>0)
      if(length(to_use_samples)<5){
        next
      }
      
      # between sample comparison (halmark pathway level variation)
      cms = base::intersect(rownames(Z_truth),rownames(Z_inferred))
      
      Z_truth_pathwayScore = singscore_result(Z_truth[cms,to_use_samples],halmark.genelist)
      Z_inferred_pathwayScore = singscore_result(Z_inferred[cms,to_use_samples],halmark.genelist)
      Total_pathwayScore =  singscore_result(simulated_bulk[cms,to_use_samples],halmark.genelist)
      
      between_sample_cor = rbind(between_sample_cor,data.frame(cell_type = paste0(ct,'->',maxCorName),
                                                               class = 'Z_inferred',
                                                               cor = pair_cor(Z_truth_pathwayScore,Z_inferred_pathwayScore,margin = 'row',cor.method = 'pearson'),
                                                               pathway = rownames(Z_truth_pathwayScore)))
      between_sample_cor = rbind(between_sample_cor,data.frame(cell_type = paste0(ct,'->',maxCorName),
                                                               class = 'Total Expr',
                                                               cor = pair_cor(Z_truth_pathwayScore,Total_pathwayScore,margin = 'row',cor.method = 'pearson'),
                                                               pathway = rownames(Z_truth_pathwayScore)))
      
      # between sample comparison (MP genes level variation)
      matched_MP_ct = MP_map[[tumorType]]$MP_ct[which(MP_map[[tumorType]]$cell_type == maxCorName)] 
      if(!is.na(matched_MP_ct)){
        Z_truth_pathwayScore = singscore_result(Z_truth[cms,to_use_samples],MP_programs[[matched_MP_ct]])
        Z_inferred_pathwayScore = singscore_result(Z_inferred[cms,to_use_samples],MP_programs[[matched_MP_ct]])
        Total_pathwayScore =  singscore_result(simulated_bulk[cms,to_use_samples],MP_programs[[matched_MP_ct]])
        
        between_sample_MP_cor = rbind(between_sample_MP_cor,data.frame(cell_type = paste0(ct,'->',maxCorName),
                                                                 class = 'Z_inferred',
                                                                 cor = pair_cor(Z_truth_pathwayScore,Z_inferred_pathwayScore,margin = 'row',cor.method = 'pearson'),
                                                                 pathway = rownames(Z_truth_pathwayScore)))
        between_sample_MP_cor = rbind(between_sample_MP_cor,data.frame(cell_type = paste0(ct,'->',maxCorName),
                                                                 class = 'Total Expr',
                                                                 cor = pair_cor(Z_truth_pathwayScore,Total_pathwayScore,margin = 'row',cor.method = 'pearson'),
                                                                 pathway = rownames(Z_truth_pathwayScore)))
        
      }
      
      
      
      # within sample comparison
      hv_genes = base::intersect(cms,hallmark_genes)
      within_sample_cor = rbind(within_sample_cor,data.frame(cell_type = paste0(ct,'->',maxCorName),
                                                             class = 'Z_inferred',
                                                             frac = simulated_frac[to_use_samples,ct],
                                                             cor = pair_cor(Z_truth[hv_genes,to_use_samples],Z_inferred[hv_genes,to_use_samples],margin = 'column',cor.method = 'spearman')) %>% set_rowNames_to_NULL())
      within_sample_cor = rbind(within_sample_cor,data.frame(cell_type = paste0(ct,'->',maxCorName),
                                                             class = 'Total Expr',
                                                             frac = simulated_frac[to_use_samples,ct],
                                                             cor = pair_cor(Z_truth[hv_genes,to_use_samples],simulated_bulk[hv_genes,to_use_samples],margin = 'column',cor.method = 'spearman')) %>% set_rowNames_to_NULL())
      
      
    }
    
    
    Z_evalu_obj[[bulk]]$between_sample_cor = between_sample_cor
    Z_evalu_obj[[bulk]]$within_sample_cor = within_sample_cor
    Z_evalu_obj[[bulk]]$between_sample_MP_cor = between_sample_MP_cor
    
    
    print(paste('**************************** finish',bulk,'***********************************'))
    

  }
  return(Z_evalu_obj)
  
}



Z_plot = function(x,fig_width = 12,fig_height = 18,plot_file_name = 'Z_performance.pdf'){
  a = x$between_sample_cor
  b = x$between_sample_MP_cor
  c = x$within_sample_cor
  
  p1 = ggplot(c,aes(frac,cor))+
    geom_point(aes(color = class),size = 1)+
    facet_wrap(~cell_type,nrow =1)+
    theme_bw()+
    ylab('spearman cor')+
    ggtitle('Halmark genes: within sample comparison')
  
  p2 = ggplot(b,aes(class,cor))+
    geom_boxplot(aes(color = class))+
    facet_wrap(~cell_type,nrow = 1)+
    theme_bw()+
    ggtitle('MP program Scores: between sample comparison')
  
  p3 = ggplot(a,aes(class,cor))+
    geom_boxplot(aes(color = class))+
    facet_wrap(~cell_type,nrow = 1)+
    theme_bw()+
    ggtitle('Halmark pathway Scores: between sample comparison')
  
  p_out=gridExtra::grid.arrange(grobs = lapply(
    list(p1,p2,p3),
    egg::set_panel_size,
    width = unit(4.5, "cm"),
    height = unit(3.75, "cm")),ncol = 1,just = "left"
  )
  
  ggsave(plot_file_name,width = fig_width,height = fig_height,plot = p_out)
}


