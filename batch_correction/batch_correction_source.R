
run_ComBat <- function(mat, batch_labels,model_matrix=NULL){
  require(sva)
  
  
  #mat <- data$df_otu_corrected
  #range(mat)
  # make continuous
  combat_predata = mat #log(mat + 1) #mat #
  input = ComBat( dat=combat_predata, batch = batch_labels,mod = model_matrix)
  return(input)
}


run_ComBat_mle <- function(mat,batch_labels){
  source(paste0(script_folder,"CBX/combatx_helper.R"))
  source(paste0(script_folder,"CBX/combat_mle.R"))
  input = ComBat_mle( dat = mat , batch = batch_labels,estimation_method="MLE")
}

run_percentile_norm <- function(mat,data,case_class, control_class){
  source(paste0(script_folder,"percentile_norm.R"))
  pernorm = percentile_norm(mat,df_meta = data$df_meta,replace_zeroes=TRUE,case_class = case_class, control_class=control_class)
  return(pernorm)
}

#' @param mat matrix you want to batch correct
#' @param data object containing $df_meta 
#' @param ref_study Reference study to align remaining samples to
run_slope_correction <- function(mat,data,ref_study){
  source(paste0(script_folder,"slope_correction.R"))
  slope_data= slope_correction_pooled(mat,data$df_meta,ref_study = ref_study)
  return(slope_data$df_otu)
}
run_limma <- function(mat,batch_labels,batch_labels2 = NULL){
  require(limma)
  input = removeBatchEffect( x=mat , batch= batch_labels,batch2 = batch_labels2)
  return(input)
}
run_bmc <- function(mat,batch_labels){
  
  
  require(dplyr)
  corrected_mat = mat
  unique_batches= unique(batch_labels)
  for( b in 1:length(unique_batches)){
    samples = colnames(mat)[batch_labels == unique_batches[b]]
    batch_mat = mat[,samples]
    corrected_mat[,samples] = sweep(mat[,samples],MARGIN = 1, rowMeans(batch_mat))
  }
  
  return(corrected_mat)
  
}

# run_dim_red_combat <- function(){
#   require(compositions)
#   ilr()
# }

pca_method <- function(input,clr_transform = FALSE,center_scale_transform =TRUE,num_pcs){
  #input = otu_data$df_otu_rel_ab
  orig_input = input
  require(compositions)
  require("bigstatsr")
  if(clr_transform){
    input = t(clr(t(input)))
  }
  if(center_scale_transform){
    input = t(scale(t(input)))
  }
  #dim(input)
  #dim(orig_input)
  
  myFBM = as_FBM(t(input), type = c("double"))
  
  
  t1 = Sys.time()
  
  #svd_result = big_SVD(myFBM,k=20,fun.scaling = big_scale())
  
  svd_result = big_SVD(myFBM,k=(num_pcs+10))

  #?big_SVD
  print(Sys.time()-t1)
  
  pca_score <-  svd_result$u %*% diag(svd_result$d)
  
  row.names(pca_score) = colnames(orig_input)
  return(list(svd_result = svd_result,pca_score=pca_score,transformed_data = input))
}

regress_out <- function(pc_scores,data,pc_index){
  
  model_residuals<-lm(as.matrix(data) ~ pc_scores[,pc_index] ) 
  extracted_residuals <- residuals(model_residuals)
  return(t(extracted_residuals))
  
}

# run_smart_sva <- function(mat, batch_labels){
#   mat = input_abundance_table
#   
#   require(SmartSVA)
#   ?EstDimRMT
#   
#   #mat = input_abundance_table
#   
#   mat_scale = t(scale(t(mat)))
#   
#   mat = mat_scale
#   
#   require(SmartSVA)
#   
#   t1 = Sys.time()
#   Y.r <- t(resid(lm(t(mat) ~ batch_labels)))
#   n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
#   print(n.sv)
#   t2= Sys.time()
# 
#   mod <- model.matrix( ~ batch_labels)
#   sv.obj <- smartsva.cpp(input_abundance_table, mod, mod0=NULL, n.sv=n.sv)
#   t2= Sys.time()
#   print(t3 -t2)
#   out_mat <- t(sv.obj$sv)
#   #row.names(out_mat) = row.names(mat)
#   colnames(out_mat) = colnames(mat)
#   
#   return(out_mat)
#   
#   out_mat[1:4,1:4]
#   dim(out_mat)
#   dim(out_mat_no_scaling)
#   out_mat_no_scaling[1:4,1:4]
#   
#   out_mat_no_scaling = out_mat
# }


#' @param mat matrix you want to batch correct
#' @param data object containing $df_meta 
#' 
#' 
# 
# mat = input_abundance_table_clr_scale
# metadata_mod= total_metadata_mod_interest
# bio_signal_formula = bio_signal_formula_interest
# num_pcs=100
# mod_ <- model.matrix( object = bio_signal_formula, data =  metadata_mod)
# mat1 = mat + 1
# sv.obj <- smartsva.cpp(dat = mat, mod = mod_, alpha = .25,
#                        mod0=NULL, n.sv=100, B = 1000, VERBOSE = T)
# 
# 
# 
# BMI=read.table(file = "~/Downloads/bmi_corrected.txt")
# library(data.table)
# kmers=as.matrix(data.frame(fread(file = "~/Downloads/kmer_table_clr_scaled.txt", header = T), row.names = 1))
# mod <- model.matrix( ~ bmi_corrected, BMI)
# 
# sv.obj <- smartsva.cpp(dat = kmers, mod = mod, alpha = .25,
#                        mod0=NULL, n.sv=100, B = 1000, VERBOSE = T)
# dim(kmers)
# 
# kmers[1:4,1:4]
# mat[1:4,1:4]



# mat = input_abundance_table
# metadata_mod = total_metadata_mod_interest
# bio_signal_formula =  bio_signal_formula_interest

run_sva <- function(mat,metadata_mod=NULL,bio_signal_formula=NULL,num_pcs = NULL){
  message("about to load smartsva")
  require(SmartSVA)
  message("finish load smartsva")
  
  
  #mat = input_abundance_table
  #metadata_mod=total_metadata_mod
  
  #mat = mat[rowVars(mat)!=0,]
  mat_scaled = mat

  
  
  message(dim(mat_scaled))
  message(dim(metadata_mod))

  message(num_pcs)
  
  if(!is.null(num_pcs)){
    n.sv = num_pcs
  }else{
    if(!is.null(metadata_mod)){
    
      
      #bio_signal_formula <- paste0("~",paste(colnames(metadata_mod), collapse = "+"))
      bio_signal_formula_resid = as.formula(paste0("t( mat_scaled)", paste0(as.character(bio_signal_formula),collapse = '')))
      
      dim(mat_scaled)
      dim(metadata_mod)
      #?smartsva.cpp
      #?smartsva.cpp
      
      #Determine number of SVs
      message("about to resid")
      Y.r <- t(resid(lm(bio_signal_formula_resid,data = metadata_mod)))
  
      message("estimating RT")
      t1 = Sys.time()
      n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1 # Very important: Add one extra dimension to compensate potential loss of 1 degree of freedom in confounded scenarios !!!
      t2= Sys.time()
      message(t2 -t1)
    }else{
      Y.r = mat_scaled
      t1 = Sys.time()
      n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1 # Very important: Add one extra dimension to compensate potential loss of 1 degree of freedom in confounded scenarios !!!
      t2= Sys.time()
      message(t2 -t1)
    }
    print("n.sv")
    print(n.sv)
  }
  

  
  # Run SVA
  detach("package:compositions", unload=TRUE)
  mod <- model.matrix( object = bio_signal_formula, data =  data.frame(metadata_mod))

  #mod <- model.matrix( object = ~ DiseaseState, data = data$df_meta)
  sv.obj <- smartsva.cpp(dat =  as.matrix(mat_scaled), mod = mod, mod0=NULL, n.sv=n.sv, B = 1000, VERBOSE = T) 
  t3= Sys.time()
  #message(t3 -t2)
  #dim(sv.obj$sv)
  #dim(mat_scaled)
  #To get corrected data run: 
  mat_scaled_corrected<- t(resid(lm(t(mat_scaled) ~ ., data=data.frame(sv.obj$sv))))
  return( list(corrected_data = mat_scaled_corrected, sv.obj=sv.obj,n.sv=n.sv))
}
