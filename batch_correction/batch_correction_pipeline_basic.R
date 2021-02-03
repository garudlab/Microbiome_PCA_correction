

args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)
#table(total_metadata$abdominal_obesity_idf_v2.y)
#table(total_metadata$diabetes_self_v2)
#table(total_metadata$diabetes_lab_v2.x)
# args = c("kmer", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "CRC", "ComBat",10,"study",1,1,"bin_crc_normal",0,"none",0,0,0,1,1)
# args = c("otu", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "CRC", "ComBat",10,"study",1,1,"bin_crc_normal",0,"none",0,0,0,1,1)

#          "CRC_thomas", "ComBat_with_biocovariates_with_seqbatch",-1,"dataset_name",1,1,"bin_crc_normal",0,"logscale",0,0,0)
# args = c("kmer", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "Thomas", "minerva",-1,"dataset_name",1,1,"bin_crc_normal",0,"none",0,0,0)

# args = c("kmer", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "AGP_max", "minerva",20,"Instrument",1,1,"bin_antibiotic_last_year",0,"clr_scale",0,0,0,1,1)
# args = c("kmer", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "AGP_max", "minerva","1","Instrument",1,1,"bin_antibiotic_last_year",0,"clr_scale",0,0,0,1,"Yes")

# args = c("otu", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "AGP_complete", "raw",10,"Instrument",1,1,"bin_antibiotic_last_year",0,"none",0,0,0,1,1)

# args = c("kmer", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "AGP_max", "PhenoCorrect",20,"Instrument",1,1,"bin_antibiotic_last_year",0,"none",0,0,0,1, "Yes")
# args = c("kmer", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "T2D", "smartsva",10,"Instrument",1,1,"bmigrp_c4_v2.x",0,"none",0,0,0,3,1)#3,1)#"1","4")
# table(total_metadata$antibiotic)

# args = c("kmer", 4, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "Hispanic", "smartsva",10,"Instrument",1,1,"mets_idf3_v2",0,"none",1,"1")
# args = c("kmer", 4, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "Hispanic", "raw",10,"extraction_robot..exp.",1,1,"bmi_v2",0,"clr_scale")

# args = c("kmer", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "PTBmeta", "refactor",10,"study",1,1,"preg_outcome",0,"clr_scale",1,"preterm")
# 
# 
# count = 1
# for(c in 1:length(colnames(total_metadata))){
#   data_na_included = as.character(total_metadata[,colnames(total_metadata)[c]])
#   data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other" 
#                    | data_na_included == '' | data_na_included == 'not applicable' | data_na_included == 'not provided'
#                    | data_na_included == 'Not applicable'| data_na_included == 'Unspecified'] = NA
#   
#   if(length(table(data_na_included))< 15 &length(table(data_na_included))> 1){
#     print(c)
#     print(colnames(total_metadata)[c])
#     print(table(data_na_included))
#     count = count + 1
#   }
# }

#table(total_metadata$income_c5_v2.x)
#table(total_metadata$income_v2.x)
#table(total_metadata$diabetes3_v2)

# args = c("otu", 6, "/u/home/b/briscoel/project-halperin/MicroBatch", "AGP_Hfilter",
#          "smartsva_clr",10,"Instrument",1, "bmi_corrected",0)

# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
num_pcs = as.integer(args[6])#5
batch_column = args[7]
save_PC_scores = as.logical(as.integer(args[8]))#TRUE
filter_low_counts = as.logical(as.integer(args[9]))
covariate_interest = args[10]
subsample_bool = 0
use_RMT = as.logical(as.integer(args[11]))
transformation = args[12]
if(length(args)> 12){
  subsample_bool = as.logical(as.integer(args[13]))
  subsample_prop =as.numeric(args[14])
  subsample_seed = as.integer(args[15])
}
if(length(args)> 15){
  label_pos_or_neg = as.integer(args[16])
  target_label = args[17]
}

# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(compositions)

script_folder = paste0(microbatch_folder,'/data_processing')
batch_script_folder = paste0(microbatch_folder, '/batch_correction')
plot_dir =paste0(microbatch_folder,'/plots/',study_name,'_k',kmer_len)
source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
#regress_out


# ============================================================================== #
# define folders
otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)
if(grepl("kmer",data_type)){
  
  output_folder = kmer_input_folder
  
}else{
  output_folder = otu_input_folder
}
if(subsample_bool){
  output_folder = paste0(output_folder, "_subsample_",as.integer(100*subsample_prop),"_seed_",subsample_seed)
  dir.create(output_folder)
}

#otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
#otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))

#total_metadata = readRDS(paste0(otu_input_folder,"metadata.rds"))


if(grepl("clr",transformation)){
  file_type = ""
}else{
  file_type = "_norm"
}
if(data_type == "kmer"){
  #dir.create(paste0(output_folder,"/",batch_column))
  
  input_folder = kmer_input_folder
  kmer_table = readRDS(paste0(kmer_input_folder,"/kmer_table", file_type,".rds"))
  
}else{
  input_folder = otu_input_folder
  otu_table = readRDS(paste0(otu_input_folder,"/otu_table", file_type,".rds"))
}
total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))
#table(total_metadata$dataset_name,total_metadata$bin_crc_normal)
#table(total_metadata$study,total_metadata$bin_crc_normal)
#table(total_metadata$Instrument,total_metadata$bin_antibiotic_last_year)

if(grepl("reprocess",study_name)){
  collection_date=as.Date(total_metadata$collection_timestamp, format="%Y-%m-%d %H:%M")
  collection_year = as.integer(format(as.Date(collection_date, format="%m/%d/%Y"), "%Y"))
  total_metadata$collection_year = collection_year
  
}


#install.packages('SmartSVA')
if(grepl("AGP",study_name)){
  new_collection_year = total_metadata$collection_year
  new_collection_year[new_collection_year < 2010] = NA
  new_collection_year[new_collection_year > 2017] = NA
  total_metadata$collection_year = new_collection_year
}



# ============================================================================== #
# cleaning of data


input_abundance_table = get(paste0(data_type,"_table"))
sum(is.na(input_abundance_table))
dim(input_abundance_table)
intersect_samples = intersect(colnames(input_abundance_table),row.names(total_metadata))

input_abundance_table = input_abundance_table[,intersect_samples]
total_metadata = total_metadata[intersect_samples,]
sum(is.na(input_abundance_table))
dim(input_abundance_table)
###@@@@@

if(subsample_bool){
  set.seed(subsample_seed)
  subsample_samples = sample(colnames(input_abundance_table),size = as.integer(subsample_prop*ncol(input_abundance_table)))
  
  
  
  if(methods_list == "ProtectPCA" | methods_list == "ProtectPCA_compare"){
    non_sample_index = which(!(colnames(input_abundance_table) %in% subsample_samples))
    test_samples  = colnames(input_abundance_table)[non_sample_index]
    train_samples = subsample_samples
    
    
    test_abundance_table = input_abundance_table[,test_samples]
    test_metadata = total_metadata[test_samples,]
    
    saveRDS(test_metadata,paste0(output_folder,"/metadata.rds"))
    write.table(test_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
    
    #train_abundance_table = input_abundance_table[,subsample_samples]
    #train_metadata = total_metadata[subsample_samples,]
    
  }else{
    input_abundance_table = input_abundance_table[,subsample_samples]
    total_metadata = total_metadata[subsample_samples,]
    saveRDS(total_metadata,paste0(output_folder,"/metadata.rds"))
    
    write.table(total_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
    
      
  }
  
  
  
}

#####@@@


# tissue_filder

if(grepl("AGP",study_name)){
  if(data_type == "kmer"){
    tissue_samples = unlist(total_metadata %>% filter(total_metadata$body_habitat.x == "UBERON:feces") %>% select(Run))
    
  }else{
    tissue_samples = unlist(total_metadata %>% filter(total_metadata$body_habitat.x == "UBERON:feces") %>% select(sample_name))
    
  }
  tissue_samples = as.character(tissue_samples)
  #length(intersect(tissue_samples,colnames(input_abundance_table)))
  input_abundance_table = input_abundance_table[,tissue_samples]
  total_metadata = total_metadata[tissue_samples,]
}

#dim(input_abundance_table)
#dim(input_abundance_table)

batch_labels = as.integer(droplevels(as.factor(total_metadata[,batch_column])))

batch_labels_dummy = to.dummy(batch_labels,"batch")
#table(total_metadata$bin_obese)
#table(batch_labels)
#"bmc","ComBat","limma",

batch_corrected_outputs = list()

if(grepl("AGP",study_name)){
  collection_date=as.Date(total_metadata$collection_timestamp, format="%Y-%m-%d")
  total_metadata$collection_date = collection_date
  collection_days = collection_date - min(collection_date,na.rm=TRUE)
  collection_month = format(as.Date(total_metadata$collection_date, format="%m/%d/%Y"), "%Y-%m")
  
  
  collection_year = as.integer(format(as.Date(collection_date, format="%m/%d/%Y"), "%Y"))
  
  batch_labels2 = as.character(collection_year)
  
  #total_metadata_mod = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",categorical_vars ="race.x",numeric_vars = "bmi_corrected")
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",categorical_vars ="race.x")
 
}else if(grepl("Thomas",study_name) | grepl("thomas",study_name) ){
  batch_labels2 = as.character(total_metadata[,'DNA_extraction_kit'])
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,binary_vars =  covariate_interest)
 
}


if(grepl("bmi",covariate_interest)){
  
  #covariate_interest = "host_body_mass_index"
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,numeric_vars =  covariate_interest)

}else if(grepl("Abx0_6",covariate_interest)){
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,categorical_vars = "antibiotic_history",
                                                     label_pos_or_neg = 0)
  new_metadata_interpretation = sapply(total_metadata_mod_interest[,1],function(x){
    if(is.na(x)){
      return(NA)
    }else if(x == "I have not taken antibiotics in the past year." | 
             x == "Year"){
      return(0)
    }else{
      return(1)
    }
    
  })
  total_metadata_mod_interest[,1] = new_metadata_interpretation
 
}else if(grepl("Abx6_12",covariate_interest)){
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,categorical_vars = "antibiotic_history",
                                                     label_pos_or_neg = 0)
  new_metadata_interpretation = sapply(total_metadata_mod_interest[,1],function(x){
    if(is.na(x)){
      return(NA)
    }else if(x == "I have not taken antibiotics in the past year."){
      return(0)
    }else if (x == "Year" | x == "6 months"){
      return(1)
    }else if ( x == "Week" | x == "Month"){
      return(NA)
    }
    
  })
  #table(new_metadata_interpretation)
  total_metadata_mod_interest[,1] = new_metadata_interpretation

}else{
  if(length(args)> 15){
    
    if( label_pos_or_neg == 3){
      total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,categorical_vars = covariate_interest,
                                                         label_pos_or_neg = label_pos_or_neg,target_label = target_label)
      total_metadata_mod_interest = model.matrix(~factor(total_metadata_mod_interest[,1]))
      colnames(total_metadata_mod_interest) = c("Intercept",paste0("Variable",c(1:(ncol(total_metadata_mod_interest)-1))))
      
    }else{
      total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,binary_vars = covariate_interest,
                                                         label_pos_or_neg = label_pos_or_neg,target_label = target_label)
      #table(total_metadata_mod_interest[,1])
      
    }
    
  }else{
    if(grepl("gest_age",covariate_interest)){
      total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,numeric_vars =  covariate_interest)
      
    }else{
      total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,binary_vars = covariate_interest)
      
    }
    
  }
  #table(total_metadata_mod_interest)

}
bio_signal_formula <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod_interest), collapse = " + ")))


bio_signal_formula_interest <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod_interest), collapse = " + ")))

# take out 0 variance rows

print(dim(total_metadata_mod_interest))
print(dim(input_abundance_table))
input_abundance_table  =input_abundance_table[,rowSums(is.na(total_metadata_mod_interest )) == 0]

dim(input_abundance_table)
sum(is.na(input_abundance_table))


batch_labels = batch_labels[rowSums(is.na(total_metadata_mod_interest )) == 0]
batch_labels2 = batch_labels2[rowSums(is.na(total_metadata_mod_interest )) == 0]


total_metadata_mod_interest= total_metadata_mod_interest[rowSums(is.na(total_metadata_mod_interest)) == 0,,drop=FALSE]


# filter out low count features
dim(input_abundance_table)
sum(is.na(input_abundance_table))

if(filter_low_counts){
  if(grepl("clr",transformation)){
    filter_at_least_two_samples_per_feature = (rowSums(input_abundance_table  > 5 ) > 2)
  }else{
    filter_at_least_two_samples_per_feature = (rowSums(input_abundance_table  > 0 ) > 2)
  }
  sum(filter_at_least_two_samples_per_feature)
  
  input_abundance_table = input_abundance_table[filter_at_least_two_samples_per_feature,]
  
}
#sum(filter_at_least_two_samples_per_feature)
#table_before_log = input_abundance_table

#input_abundance_table = table_before_log
#table_before_log[1:5,1:5]
#rowSums(table_before_log)[1:5]
if(grepl("log",transformation)){
  print("LOGTIME")
  input_abundance_table = log(input_abundance_table+1)
  
}
dim(input_abundance_table)
sum(is.na(input_abundance_table))

if(grepl("clr",transformation)){
  input_abundance_table = t(clr(t(input_abundance_table)))
  input_abundance_table = data.frame(input_abundance_table)
  input_abundance_table = as.matrix(input_abundance_table)
  
}


if(grepl("scale",transformation)){
  input_abundance_table= t(scale_custom(t(input_abundance_table)))
  
}

input_abundance_table = input_abundance_table[rowVars(as.matrix(input_abundance_table)) > 10e-10 ,]


if(!grepl("reprocess",study_name) & grepl("AGP",study_name)){
  total_metadata_mod2 = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",
                                             categorical_vars =c("bin_omnivore_diet","bin_antibiotic_last_year"))
  bio_signal_formula2 <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod2), collapse = " + ")))
  
  
  
}

print("dimensions")
print(dim(input_abundance_table))




for(m in 1:length(methods_list)){
  sv_object_output = c()
  
  batch_corrected_output  = c()
  
  input_abundance_table_mod = c()
  batch_corrected_output1 = c()
  batch_corrected_output1_mod = c()
  batch_labels2_mod = c()
  
  if(methods_list[m] == "bmc"){
    batch_corrected_output = run_bmc(mat = input_abundance_table, batch_labels)
    #length(batch_labels)
    #dim(input_abundance_table)
  }else if(methods_list[m] == "raw"){
    
    batch_corrected_output = input_abundance_table
    pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = TRUE,num_pcs = 10 )
    sv_object_output =  pca_res
    # if(save_PC_scores){
    #   sv_object_output =  pca_res
    #   row.names(sv_object_output$sv) = colnames(input_abundance_table)
    #   
    #   saveRDS(sv_object_output, paste0(output_folder ,"/protect_",covariate_interest,"/SVs_",methods_list[m],extra_file_name,".rds"))
    #   
    #   # write.table(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder,"/",batch_column,"/BatchCorrected_",methods_list[m],extra_file_name,".txt"),
    #   #             sep = "\t",quote = FALSE)
    #   
    # }
    
  }else if(methods_list[m] == "clr"){
    require(compositions)
    batch_corrected_output = input_abundance_table_clr
    
  }else if(methods_list[m] == "ilr"){
    require(compositions)
    batch_corrected_output = t(ilr(t(input_abundance_table)))
    
  }else if(methods_list[m] == "ComBat"){
   
    batch_corrected_output = run_ComBat(mat = input_abundance_table, batch_labels)
    
    
  }else if(methods_list[m] == "ComBatLog"){
    batch_corrected_output = run_ComBat(mat = log(input_abundance_table), batch_labels)
    
  }else if(methods_list[m] == "ComBatLog_with_batch2"){
    batch_corrected_output1 = run_ComBat(mat = log(input_abundance_table), batch_labels)
    
    batch_corrected_output1_mod =  batch_corrected_output1[,!is.na(batch_labels2)]
    batch_labels2_mod = batch_labels2[!is.na(batch_labels2)]
    batch_corrected_output = run_ComBat(mat = batch_corrected_output1_mod, batch_labels = batch_labels2_mod)
    
  }else if(methods_list[m] == "ComBat_with_batch2"){
    batch_corrected_output1 = run_ComBat(mat = input_abundance_table, batch_labels)
    
    batch_corrected_output1_mod =  batch_corrected_output1[,!is.na(batch_labels2)]
    batch_labels2_mod = batch_labels2[!is.na(batch_labels2)]
    batch_corrected_output = run_ComBat(mat = batch_corrected_output1_mod, batch_labels = batch_labels2_mod)

    
  }else if(methods_list[m] == "ComBat_with_biocovariates"){
    total_metadata_mod = total_metadata_mod_interest
    input_abundance_table_mod = input_abundance_table[,rowSums(is.na(total_metadata_mod)) == 0]
    metadata_mod= total_metadata_mod[rowSums(is.na(total_metadata_mod)) == 0,]
    batch_labels_mod = batch_labels[rowSums(is.na(total_metadata_mod)) == 0]
    mod <- model.matrix( bio_signal_formula, data = total_metadata_mod)
    batch_corrected_output = run_ComBat(mat = input_abundance_table_mod, batch_labels = batch_labels_mod,mod = mod)
    
    
  }else if(methods_list[m] == "ComBat_with_biocovariates_with_seqbatch"){
    #total_metadata_mod1 = total_metadata_mod
    #total_metadata_mod1[total_metadata_mod1 == "African American"] = NA
    #metadata_mod$race.x= as.factor(as.character(metadata_mod$race.x))
    
    total_metadata_mod = total_metadata_mod_interest
    input_abundance_table_mod = input_abundance_table[,rowSums(is.na(total_metadata_mod)) == 0]
    metadata_mod= total_metadata_mod[rowSums(is.na(total_metadata_mod)) == 0,,drop=FALSE]
    batch_labels_mod = batch_labels[rowSums(is.na(total_metadata_mod)) == 0]
    batch_labels2_mod = batch_labels2[rowSums(is.na(total_metadata_mod)) == 0]
    head(total_metadata_mod)
    mod <- model.matrix( bio_signal_formula, data = metadata_mod)
    
    batch_corrected_output1 = run_ComBat(mat = input_abundance_table_mod, batch_labels = batch_labels_mod,mod = mod)
    
    tab_batch2 = table(batch_labels2_mod)
    lonely_batches = names(tab_batch2[tab_batch2 == 1])
    batch_corrected_output1 = batch_corrected_output1[,!(batch_labels2_mod %in% lonely_batches) & !is.na(batch_labels2_mod)]
    metadata_mod= metadata_mod[!(batch_labels2_mod %in% lonely_batches)& !is.na(batch_labels2_mod),,drop=FALSE]
    batch_labels2_mod = batch_labels2_mod[!(batch_labels2_mod %in% lonely_batches)& !is.na(batch_labels2_mod)]
    mod <- model.matrix( bio_signal_formula, data = metadata_mod)
    
    batch_corrected_output2 = run_ComBat(mat = batch_corrected_output1, batch_labels2_mod,mod = mod)
    batch_corrected_output = batch_corrected_output2
    
  }else if(methods_list[m] == "minerva"){
    set.seed(0)
    pca_res = 0
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva_clr",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
    }else{
      num_factors = num_pcs
    }
    
    if(num_pcs == -1){
      
      pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE,num_pcs = 20 )
      
      
      sv_object_output = pca_res
      pca_data = pca_res$pca_score
      transformed_data = t(pca_res$transformed_data)
      
      
      
      #Randomly shuffle the data
      shuffle_samples <-sample(row.names(transformed_data))
      
      #Create 10 equally size folds
      n_cv = 3
      original_time = Sys.time()
      folds <- cut(seq(1,length(shuffle_samples)),breaks=n_cv,labels=FALSE)
      
      calibration_list_cv = list()
      
      #Perform 10 fold cross validation
      
      for(i in 2:n_cv){
        calibration_list_cv[[i]] = list()
        #Segement your data by fold using the which() function 
        testIndexes <- which(folds==i,arr.ind=TRUE)
        testSamples = shuffle_samples[testIndexes]
        trainSamples = shuffle_samples[-testIndexes]
        test_pca_data <-  pca_data[testSamples, ]
        train_pca_data <- pca_data[trainSamples, ]
        
        test_transformed_data <- transformed_data[testSamples,]
        train_transformed_data<- transformed_data[trainSamples,]
        
        
        
        calibration_list = list()
        for(num_factors_it in 1:20){
          calibration_list[[num_factors_it]] = regress_out(train_pca_data,data=train_transformed_data,pc_index = c(1:num_factors_it))
          
          
        }
        
        library(randomForest)
        # Perform training:
        max_accuracy = c()
        for(num_factors_it in 1:20){
          #num_factors_it = 1
          corrected_data = calibration_list[[num_factors_it]]
          phenotype = total_metadata_mod_interest[trainSamples,]
          not_na_samples = !is.na(phenotype)
          corrected_data = corrected_data[,not_na_samples]
          phenotype = phenotype[not_na_samples]
          
          
          # time1 = Sys.time()
          # rf_classifier = randomForest(phenotype ~ ., data=t(corrected_data) , ntree=100, mtry=2, importance=TRUE)
          # print(Sys.time() - time1)
          # predict(rf_classifier, t(corrected_data))
          # 
          # library(e1071)
          # time1 = Sys.time()
          # nb_classifier = naiveBayes(phenotype ~ ., data=t(corrected_data) )
          # print(Sys.time() - time1)
          # predict(nb_classifier,t(corrected_data))
          
          
          require(caret)
          # set up 10-fold cross validation procedure
          train_control <- trainControl(
            method = "cv", 
            number = 3
          )
          
          time1 = Sys.time()
          # train model
          nb.m1 <- train(
            x = t(corrected_data),
            y = phenotype,
            method = "nb",
            trControl = train_control
          )
          print(Sys.time() - time1)
          
          max_accuracy = c(max_accuracy,max(nb.m1$results$Accuracy))
          
          
        }
        
        print(length(max_accuracy))
        opt_num_pcs = which.max(max_accuracy)
        
        calibration_list_cv[[i]][["all_train_performance"]] = max_accuracy
        calibration_list_cv[[i]][["max_train_performance"]]  = max(max_accuracy)
        calibration_list_cv[[i]][["optimal_PC"]]  = which.max(max_accuracy)
        
        
        test_corrected = regress_out(test_pca_data,data=test_transformed_data,pc_index = c(1:opt_num_pcs))
        phenotype_test = total_metadata_mod_interest[testSamples,]
        not_na_samples = !is.na(phenotype_test)
        test_corrected  = test_corrected [,not_na_samples]
        phenotype_test = phenotype_test[not_na_samples]
        
        require(klaR)
        #full_data = 
        
        nb_test_pred = predict(nb.m1$finalModel,newdata = data.frame(t(test_corrected)))
        #length(nb_test_pred$class)
        #length(phenotype_test)
        #dim(test_corrected)
        #dim(corrected_data)
        #table(phenotype_test)
        #table(total_metadata$bin_crc_normal)
        conf_mat = confusionMatrix(nb_test_pred$class, phenotype_test)
        calibration_list_cv[[i]][["test_performance"]] = conf_mat$overall["Accuracy"]
        #Use the test and train data partitions however you desire...
      }
      print("Whole time")
      print(Sys.time() - original_time)
      do.call(rbind,calibration_list_cv)
      
      
      
      ####
      
      
      
      
      
    }else{
      pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE,num_pcs = num_factors )
      sv_object_output = pca_res
      
      
      batch_corrected_output = regress_out(pca_res$pca_score,data=t(pca_res$transformed_data),pc_index = c(1:num_factors))
      
    }
    # pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE,num_pcs = num_factors )
    # 
    # sv_object_output = pca_res
    # 
    # 
    # batch_corrected_output = regress_out(pca_res$pca_score,data=t(pca_res$transformed_data),pc_index = c(1:num_factors))
    # 
  }else if(methods_list[m] == "ProtectPCA"){
    
    
    set.seed(0)
    pca_res = 0
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva_clr",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
    }else{
      num_factors = num_pcs
    }
    
    
    pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE,num_pcs = num_factors )
    
    sv_object_output = pca_res
    
    all_samples = colnames(input_abundance_table)
    
    
    train_samples_final = intersect(train_samples,all_samples)
    test_samples_final = intersect(test_samples,all_samples)
    
    
    train_samples_index = which(colnames(input_abundance_table) %in% train_samples_final)
    train_metadata_mod_interest = total_metadata_mod_interest[train_samples_index,,drop=FALSE]
    train_pca_score = pca_res$pca_score[train_samples_final,]
    
    
    test_samples_index = which(colnames(input_abundance_table) %in% test_samples_final)
    test_pca_score = pca_res$pca_score[test_samples_final,]
    
    #$#$#$
    
    if(grepl("bmi",covariate_interest)){
      corr_minerva = cor(train_metadata_mod_interest,train_pca_score[,1:num_factors])
      protected_pc = which.max(abs(corr_minerva))
      
    }else{
      corr_minerva = cor(as.integer(train_metadata_mod_interest[,1])-1,train_pca_score[,1:num_factors])
      protected_pc = which.max(abs(corr_minerva))
    }
    
    #$#$#
    
    batch_corrected_output = regress_out(test_pca_score,data=t(pca_res$transformed_data[,test_samples_final]),pc_index =c(1:(protected_pc-1)))
    
  }else if(methods_list[m] == "ProtectPCA_compare"){
    set.seed(0)
    pca_res = 0
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva_clr",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
    }else{
      num_factors = num_pcs
    }
    
    
    pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE,num_pcs = num_factors )
    
    sv_object_output = pca_res
    
    all_samples = colnames(input_abundance_table)
    
    
    train_samples_final = intersect(train_samples,all_samples)
    test_samples_final = intersect(test_samples,all_samples)
    
    
    train_samples_index = which(colnames(input_abundance_table) %in% train_samples_final)
    train_metadata_mod_interest = total_metadata_mod_interest[train_samples_index,,drop=FALSE]
    train_pca_score = pca_res$pca_score[train_samples_final,]
    
    
    test_samples_index = which(colnames(input_abundance_table) %in% test_samples_final)
    test_pca_score = pca_res$pca_score[test_samples_final,]
    
    
    #$#$#
    
    batch_corrected_output = regress_out(test_pca_score,data=t(pca_res$transformed_data[,test_samples_final]),pc_index =c(1:num_factors))
    
  }else if(methods_list[m] == "limma"){
    #table(batch_labels2)
    batch_corrected_output = run_limma(mat = input_abundance_table, batch_labels)
    #batch_corrected_outputs[["limma"]] = batch_corrected_output
  }else if(methods_list[m] == "limma_batch2"){
    
    #dim(total_metadata)
    #dim(input_abundance_table)
    #length(batch_labels)
    #length(collection_days)
    
    input_abundance_table_mod =  input_abundance_table[,!is.na(batch_labels2)]
    batch_labels_mod =  batch_labels[!is.na(batch_labels2)]
    batch_labels2_mod = batch_labels2[!is.na(batch_labels2)]
    
    batch_corrected_output = run_limma(mat = input_abundance_table_mod, batch_labels = batch_labels_mod,batch_labels2 = batch_labels2_mod)
    #input = removeBatchEffect( x=input_abundance_table , batch= batch_labels,batch2 = collection_days,covariates = )
    #cbind()
    
  }else if(methods_list[m] == "smartsva"){
    print("about to start smartsva")
    if(use_RMT){
      
      sva_result= run_sva(mat = input_abundance_table, metadata_mod=total_metadata_mod_interest,bio_signal_formula = bio_signal_formula_interest,num_pcs=NULL)
    }else{
      sva_result= run_sva(mat = input_abundance_table, metadata_mod=total_metadata_mod_interest,bio_signal_formula = bio_signal_formula_interest,num_pcs=num_pcs)
      
      
      
    }
    #sva_result$corrected_data[1:4,1:4]
    #total_metadata_mod_interest
    
    
    #test =readRDS("~/Downloads/sva_numbers.rds")
    #all(test$corrected_data == sva_result$corrected_data)
    
    #table(as.factor(total_metadata_mod_interest))
    svobj = sva_result$sv.obj
    sv_object_output =  svobj
    row.names(sv_object_output$sv) = colnames(input_abundance_table)
    
    
    batch_corrected_output = sva_result$corrected_data
    
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_",methods_list[m],".txt"))
      writeLines(as.character(sva_result$n.sv), fileConn)
      message(paste0("NUM sv ", sva_result$n.sv))
    }
    num_factors = as.integer(sva_result$n.sv)
    
  }else if(methods_list[m ] == "refactor"){
    
    # rs = rowSums(input_abundance_table)
    # rv = rowVars(input_abundance_table)
    # sum(rv < 10e-9)
    # input_abundance_table_rf = input_abundance_table[(rv > 10e-4),]
    # dim(input_abudance_table_rf)
    # dim(input_abundance_table)
    # 
    
    require(TCA)
    
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
      
    }else{
      num_factors = num_pcs
    }
    
    
    sd_list = sqrt(rowSds(input_abundance_table))
    print(paste0("mean(sd):",mean(sd_list)))
    #hist(sd_list,breaks = 100)
    if(grepl(study_name,"AGP")){
      refactor_res = refactor(input_abundance_table, k=num_factors,sd_threshold = 0.001)
    }else{
      refactor_res = refactor(input_abundance_table, k=num_factors,sd_threshold = 0.02)
    }
    
    #write.table(input_abundance_table,"~/Downloads/RefactorExample.txt",quote = FALSE,sep = "\t")
    
    #sum(rowSums(input_abundance_table[,1:30])==0)
    
    RC = refactor_res$scores
    mat_scaled_corrected<- t(resid(lm(t(input_abundance_table) ~ ., data=data.frame(RC))))
    
    #row.names(refactor_res$scores) = colnames(input_abundance_table)
    sv_object_output= refactor_res
    
    #row.names(sv_object_output$scores) = row.names(input_abundance_table)
    batch_corrected_output = mat_scaled_corrected
  }else if(methods_list[m ] == "refactor_protect"){
    
    require(TCA)
    
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
      
    }else{
      num_factors = num_pcs
    }
    
    #dim(total_metadata_mod_interest)
    
    
    refactor_res = refactor(input_abundance_table, k=num_factors,C=total_metadata_mod_interest, C.remove =TRUE)
    
    RC = refactor_res$scores
    mat_scaled_corrected<- t(resid(lm(t(input_abundance_table) ~ ., data=data.frame(RC))))
    row.names(refactor_res$scores) = colnames(input_abundance_table)
    
    sv_object_output= refactor_res
    #row.names(sv_object_output$scores) = row.names(input_abundance_table)
    batch_corrected_output = mat_scaled_corrected
  }else if(methods_list[m ] =="PhenoCorrect"){
    
    
    
    
    batch_labels_factor = factor(batch_labels)
    batch_mat = model.matrix(~batch_labels_factor )
    phen_correct<- t(resid(lm(as.numeric(total_metadata_mod_interest[,1]) ~ batch_mat)))
    phen_correct = round(phen_correct[1,],10)
    # wrong approach:
    new_metadata = total_metadata[colnames(input_abundance_table),]
    new_metadata[,covariate_interest] = phen_correct
    batch_corrected_output = input_abundance_table
    
    
    output_folder = paste0(output_folder, "_",methods_list[m])
    dir.create(output_folder)
    saveRDS(new_metadata,paste0(output_folder,"/metadata.rds"))
    write.table(new_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
    
    #phen_correct_matrix = as.matrix(phen_correct)  
    #colnames(phen_correct_matrix) = colnames(input_abundance_table)
    
    #batch_corrected_output = phen_correct_matrix
    
  }else if(methods_list[m ] == "DataAugmentation"){
    batch_labels_factor = factor(batch_labels)
    batch_mat = model.matrix(~batch_labels_factor )
    input_abundance_table_aug = rbind(input_abundance_table,t(batch_mat[,2:ncol(batch_mat)]))
    batch_corrected_output = input_abundance_table_aug
    
    output_folder = paste0(output_folder, "_",methods_list[m])
    dir.create(output_folder)
    saveRDS(total_metadata,paste0(output_folder,"/metadata.rds"))
    write.table(total_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
    
  }else if(methods_list[m ] =="DomainCorrect"){
    
    batch_labels_factor = factor(batch_labels)
    batch_mat = model.matrix(~batch_labels_factor )
    message("Dim batch mat")
    message(dim(batch_mat))
    
    input_table_correct<- t(resid(lm(t(input_abundance_table) ~ batch_mat)))
    
    #output_folder = paste0(output_folder, "_",methods_list[m])
    #dir.create(output_folder)
    #saveRDS(new_metadata,paste0(output_folder,"/metadata.rds"))
    #write.table(new_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
    
    batch_corrected_output = input_table_correct 
    
    
  }else if(methods_list[m ] =="PredDomainPheno"){
    
    domain_pheno = paste0(batch_labels,"_",total_metadata_mod_interest[,1])
    
    new_metadata = total_metadata[colnames(input_abundance_table),]
    new_metadata$domain_pheno = domain_pheno 
    batch_corrected_output = input_abundance_table
    
    
    output_folder = paste0(output_folder, "_",methods_list[m])
    dir.create(output_folder)
    saveRDS(new_metadata,paste0(output_folder,"/metadata.rds"))
    write.table(new_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
    
    #phen_correct_matrix = as.matrix(phen_correct)  
    #colnames(phen_correct_matrix) = colnames(input_abundance_table)
    
    #batch_corrected_output = phen_correct_matrix
    
  }
  
  if(grepl("Abx",covariate_interest)){
    
    output_folder = paste0(output_folder, "_",covariate_interest)
    dir.create(output_folder)
    
    new_metadata = total_metadata[colnames(input_abundance_table),]
    old_colnames = colnames(new_metadata)
    new_metadata = cbind(new_metadata,total_metadata_mod_interest[,1])
    colnames(new_metadata) = c(old_colnames,covariate_interest)
    batch_corrected_output = input_abundance_table
    saveRDS(new_metadata,paste0(output_folder,"/metadata.rds"))
    write.table(new_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
  }
  #names(batch_corrected_outputs)
  batch_corrected_outputs[[methods_list[m]]] =  batch_corrected_output
  #names(batch_corrected_outputs)
  #batch_corrected_outputs[["smartsva_no_scale"]] = out_mat_no_scaling
  #batch_corrected_outputs[["smartsva_scale"]] = out_mat
  
  # make file_name
  dir.create(paste0(output_folder,"/protect_",covariate_interest))
  
  
  extra_file_name= ""
  if(grepl("pca",methods_list[m]) |grepl("refactor",methods_list[m]) |grepl("sva",methods_list[m]) | grepl("minerva",methods_list[m]) |grepl("Protect",methods_list[m]) ){
    extra_file_name = paste0(extra_file_name,"_first",num_factors)
  }
  extra_file_name = paste0(extra_file_name,"filter_",filter_low_counts, "_trans_",transformation)
  
  old_col_names = colnames(batch_corrected_outputs[[methods_list[m]]])
  if(substr(old_col_names[1],1,1) == "X"){
    new_col_names = gsub("X","",old_col_names)
    temp = batch_corrected_outputs[[methods_list[m]]]
    colnames(temp) = new_col_names
    batch_corrected_outputs[[methods_list[m]]] = temp
  }
  write.table(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder,"/protect_",covariate_interest,"/BatchCorrected_",methods_list[m],extra_file_name,".txt"),
              sep = "\t",quote = FALSE)
  saveRDS(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder ,"/protect_",covariate_interest,"/BatchCorrected_",methods_list[m],extra_file_name,".rds"))
  if(save_PC_scores){
  
    saveRDS(sv_object_output, paste0(output_folder ,"/protect_",covariate_interest,"/SVs_",methods_list[m],extra_file_name,".rds"))
    
    # write.table(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder,"/",batch_column,"/BatchCorrected_",methods_list[m],extra_file_name,".txt"),
    #             sep = "\t",quote = FALSE)
    
  }
  
  
  
  
  
}

# write.table(input_abundance_table_clr_scale,paste0(output_folder,"/",batch_column,"/kmer_table_clr_scaled.txt"),
#             sep = "\t",quote = FALSE)
# write.table(total_metadata_mod_interest,paste0(output_folder,"/",batch_column,"/bmi_corrected.txt"),
#             sep = "\t",quote = FALSE)
#write.table(input_abundance_table ,paste0(kmer_input_folder ,"/BatchCorrected_raw.txt"),sep = "\t",quote = FALSE)


#smartsva = readRDS("~/Downloads/batch_correction_outputs.rds")
#smartsva$smartsva_no_scale
# 
# data$df_meta$S = as.integer(as.factor(data$df_meta$study))
# data$df_meta$DS = as.integer(as.factor(data$df_meta$DiseaseState))
# 
# library(variancePartition)
# dir.create(plot_path)
# dir.create(paste0(plot_path,"VariancePartition"))
# for(m in 1:length(methods_list)){
#   pre = batch_corrected_outputs[[m]]
#   pre = pre[rowSums(pre)!=0,]
#   varpar=fitExtractVarPartModel(formula = ~  1 +S + DS, exprObj = pre, data = data$df_meta)
#   pdf(paste0(plot_path,"VariancePartition/",methods_list[m],".pdf"))
#   plot(plotVarPart(varpar,main=methods_list[m]))
#   dev.off()
# }

# sum(is.na(data$df_otu_clr))
# colSums(data$df_otu_clr)
# pca_rel = prcomp(t(data$df_otu_rel_ab))
# pca_rel$x[,2]
# source(paste0(main_folder,"ForHoffman/plotting_source.R"))
# 
# pca_res = pca_fn(data$df_otu_rel_ab,sample_column_true=TRUE,label_strings=data$df_meta$study,
#        filename=paste0(plot_path,"PCA/","rel_ab"),title="Pca on rel ab",acomp_version = FALSE)

#length(sv_object_output$svd_result$d)
#plot(sv_object_output$svd_result$d)

