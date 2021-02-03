#!/usr/bin/env python
# coding: utf-8

# In[1]:


#args =  ["classifier.py","/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "CRC_k6&CRC_k7", "kmer", "BatchCorrected",  "bin_crc_normal",1,0,20,1,1] 

#python MINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc "CRC_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 20 1 1
# python MINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc "AGP_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 20 1 1

# In[2]:

# from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
# data = np.random.rand(100,80)
# labels = np.random.binomial(size=100, n=1, p= 0.5)

# parameter_dict = {'n_estimators':[10],'criterion': ['entropy'],'min_samples_leaf': [10],'max_features':[0.3],'min_samples_split': [5],'max_depth':[1]}
# rf = RandomForestClassifier()
# clf = GridSearchCV(rf, parameter_dict,scoring="roc_auc")
# clf.fit(data, labels)
# clf.best_params_
# y_pred = clf.predict_proba(data)
# y_true = labels
# already_trained_auc = roc_auc_score(y_true = y_true, y_score = y_pred[:,1])

# already_trained_auc
# reducing n_estimators reduced AUC on train, increasing min_samples_leaf reduced AUC on train
use_validation = True 



import sys
import pandas as pd
import utils
import numpy as np
from sklearn.preprocessing import StandardScaler, normalize, OneHotEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn import model_selection 
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split, LeaveOneGroupOut
import statsmodels.formula.api as sm
from sklearn.metrics import roc_auc_score
from collections import Counter
from timeit import default_timer as timer
import pickle



args = sys.argv
print(args)

whole_time_start = timer()
# In[3]:


greater_folder = args[1] # what folder do you save your different datasets in
study_names = args[2].split("&")  # what is the name of the dataset (tells the program which folder to check)
data_type = args[3] # type of data. kmer vs OTU

prefix_name = args[4] # what is the prefix of the file name
column_of_interest = args[5] # what is the phenotype you are predicting (use the same name in the column of the metadata you want to predict), this programs reads from metadata.txt

norm_input = bool(int(args[6]))
map_with_accession = bool(int(args[7]))
n_repeats_input = int(args[8])
num_pcs = 20
num_pcs = int(args[9])
special_name = args[10]

perform_MINERVA = int(args[11]) # 1: minerva, #2 predict with PC's directly

spec_label_scheme = bool(int(args[12]))

label_pos_or_neg = int(args[13]) # do you want to treat CRC as positive class or negative class? 
target_label = args[14] # phenotype representing positive class or negative class? eg. CRC eg. H
print(target_label)
if not spec_label_scheme:
    label_pos_or_neg = 1
    target_label = 1


n_estimators_input = int(args[15])
criterion_input = str(args[16])
min_samples_leaf_input = int(args[17])
max_features_input = float(args[18])
min_samples_split_input = int(args[19])
max_depth_input = int(args[20])
train_it_input = int(args[21])

if len(args) > 22:
    print("LODO time?")
    bool_lodo = bool(int(args[22]))
    lodo_group = args[23]
    bool_sep_pc = bool(int(args[24])) # perform PCA on training set only and apply eigen vectors to test set
else:
    bool_lodo = False




file_output_string  = "GRID_nest" + str(n_estimators_input) + "_cri" + str(criterion_input) + "_min" + str(min_samples_leaf_input) + \
    "_max" + str(max_features_input) + "_msp" + str(min_samples_split_input) + "_mad" + str(max_depth_input) + "_trainit" + str(train_it_input)


use_domain_pheno = False # for when running raw to compare to domain pheno
if data_type == "otu" or data_type == "kmer":
    output_folders = [greater_folder + "/data/" + study_name + "/" for study_name in study_names]
    data_folders = [greater_folder + "/data/" + study_name + "/" for study_name in study_names] 
    metadata_folder =   greater_folder + "/data/" + study_names[0] + "/"  
else:
    output_folders = [greater_folder + "/data/" + study_name + "/" for study_name in study_names]
    data_folders = [greater_folder + "/data/" + study_name + "/" + "protect_" + column_of_interest + "/" + prefix_name + "_"  for study_name in study_names] 
    metadata_folder = greater_folder + "/data/" + study_names[0] + "/" 



#########################################################################
###### COMMENTARY: function definitions ######
#########################################################################




def pca_regression(y,X):
    model = sm.OLS(y,X)
    results = model.fit()
    predictedValues = results.predict()
    residuals = y - predictedValues
    return(residuals)

def RF_grid_search(data,labels,param_dict):
    rf = RandomForestClassifier()
    clf = GridSearchCV(rf, param_dict,scoring="roc_auc")
    clf.fit(data, labels)

    print("Best parameters set found on development set:")
    best_params = clf.best_params_
    return clf,best_params



def RF_cv(data,labels,param_dict):
    clf = RandomForestClassifier(max_depth=5, random_state=0,n_estimators = param_dict['n_estimators'],            criterion = param_dict['criterion'],min_samples_leaf = param_dict['min_samples_leaf'],                           max_features = param_dict['max_features'])
    results = cross_val_score(clf,X=data,y=labels,scoring="roc_auc")
    return(results)


def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

#########################################################################
###### COMMENTARY: load data from your k-mer matrix, load metadata ######
#########################################################################



# Regress out PCs
all_datasets_dict =  dict()
# For each dataset (kmer size)

start = timer()


# filter out bad metadata
metadata = pd.read_csv(metadata_folder + "metadata.txt",delimiter="\t")
if "AGP" in study_names[0]:
    tissue_samples = metadata.index[metadata['body_habitat.x'] == "UBERON:feces"]
    metadata = metadata.loc[tissue_samples]

#print(Counter(metadata[column_of_interest]))
print("new metadata shape after feces filter")
print(metadata.shape)

if spec_label_scheme :
    if label_pos_or_neg == 1:
        print("positive")
        metadata[column_of_interest] = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'],pos_labels =[target_label])
    elif label_pos_or_neg == 0:
        metadata[column_of_interest] = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'],neg_labels =[target_label])

print(Counter(metadata[column_of_interest]))

non_nan_samples = metadata.index[np.invert(np.isnan(metadata[column_of_interest]))]


# Random forest stuff
n_splits = 5
n_repeats = n_repeats_input
import random
print ("Random number with seed 30")
random.seed(30)

rskf = model_selection.RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)
logo = LeaveOneGroupOut()
parameter_dict = {'n_estimators':[n_estimators_input],'criterion': [criterion_input],\
'min_samples_leaf': [min_samples_leaf_input],'max_features':[max_features_input],\
'min_samples_split': [min_samples_split_input],'max_depth':[max_depth_input]}



for d in range(len(study_names)): # range(1):#
    ### COPY INTO loop
    feature_table_dict = utils.load_feature_table([data_folders[d]],data_type = data_type)
    feature_table = feature_table_dict[0]
    

    if norm_input:
        temp = pd.DataFrame(normalize(feature_table.transpose(), axis = 1, norm = 'l1').transpose())
        temp.index = feature_table.index
        temp.columns = feature_table.columns
        feature_table = temp

    if "AGP" in study_names[0]:
        print("feature table columns")
        print(feature_table.columns)
        tissue_samples = intersection(feature_table.columns,tissue_samples)
        feature_table = feature_table[tissue_samples]



        


    print(Counter(metadata[column_of_interest]))

    print("pos label")
    print(target_label)

    #########################################################################
    ###### COMMENTARY:  efining labels and binarize if not already     ######
    #########################################################################

    non_nan_samples = intersection(feature_table.columns,non_nan_samples)
    #print(non_nan_samples)

    feature_table = feature_table[non_nan_samples]

    metadata_labels = metadata.loc[non_nan_samples]

    labels = metadata_labels[column_of_interest]



        
    # outline
    # for data_table:
        # test train split
        # for train
        # for PC number:
            # regress out PCs
            # for grid cell:
                # get accuracy:
            # get max accuracy
        # get max accuracy across PCs
        # for test
        # get the test accuracy with best parameters
        # save best grid and best PCs

    ###########################
    #####  Preparing Data #####
    ###########################

    if perform_MINERVA == 0 or perform_MINERVA == 3 or perform_MINERVA == 4: # 3 is for data augment , 4 is for correcting oTU with domains
        feature_table_np = np.array(feature_table)
        labels_np = np.array(labels)
        dataset_start= timer()
        results_dict = dict()
        X = feature_table_np.transpose()
        y = labels_np
        na_mask = pd.isna(y)
        # samples by features
        X = X[~na_mask,:] # get rid of samples with na labels
        y = y[~na_mask]
        metadata_labels_temp = metadata_labels.loc[~na_mask,:]
        # for each test train split in 5 fold cross validation
        train_it = 0
        print("before adding dummies")
        print(X[295:300,(X.shape[1]-5):X.shape[1]])
        groups = np.array(metadata_labels_temp[lodo_group])
        groups_one_hot = pd.get_dummies(groups,drop_first=True)
        if perform_MINERVA == 3:
            #X =  [[3,4,5],[7,8,9],[6,5,6],[8,8,8]]
            X = np.append(X, groups_one_hot, 1)
            print("after adding dummies")
            print(X[295:300,(X.shape[1]-5):X.shape[1]])
        if perform_MINERVA == 4:
            #X = pca_regression(X_train,pc_scores_train[:,0:p])
            X = pca_regression(X, groups_one_hot)
            print("after regressing domain")
            print(X[295:300,(X.shape[1]-5):X.shape[1]])
        if bool_lodo:
            print("lodo time")
            logo = LeaveOneGroupOut()
            splitter = logo.split(X, y, groups)

        else:
            splitter = rskf.split(X, y)
        for train_index, test_index in splitter:   
            #print("len train index")
            #print(train_index[0:5]) 
            #print(len(train_index))
            #print(len(test_index))

            if train_it == train_it_input: 
                results_dict = dict()
                results_dict['train_best_params'] = dict()
                results_dict['train_auc_trained'] = []
                results_dict['mean_train_cv_auc'] = []
                results_dict['mean_test_cv_auc'] = []
                results_dict['test_auc_trained'] = []
                results_dict['val_auc_trained']= []         
                test_train_start = timer()
                X_train, X_test = X[train_index,], X[test_index,]
                y_train, y_test = y[train_index], y[test_index]
                if use_validation:
                    n_splits = 5
                    n_repeats = 1
                    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train,test_size=0.30, random_state=1) 
                    #print("First 5 xval")
                    #print(X_val[0:5])
                
                results_dict["number samples"] = []
                results_dict["number samples"].append(X_train.shape[0])
                # for each PC we regress out 
                
                    
                # perform grid search on train
                best_train_model, best_params = RF_grid_search(X_train, y_train,parameter_dict)
                

                pickle.dump(best_train_model, open( output_folders[d] + special_name + file_output_string  + "_grid.pkl", "wb" ) )
                test_train_end = timer()
                print("Finished one test train split")
                print(test_train_end - test_train_start)
                
            train_it += 1
            



    elif perform_MINERVA == 1:
        # get PC scores
        pca = PCA(n_components=num_pcs,svd_solver='randomized')
        if not bool_lodo:
            # do the PCA thing
            temp = feature_table.transpose()
            pca.fit(temp)
            print("Time for pca cal")
            pca_start= timer()
            print("pca start clock")
            print(pca_start)
            pc_table = pca.transform(temp)  
            print(timer() - pca_start )
        feature_table_np = np.array(feature_table)
        labels_np = np.array(labels)
         
        dataset_start= timer()
        results_dict = dict()
        X = feature_table_np.transpose()
        y = labels_np
        
        na_mask = pd.isna(y)
        print("Xshape before")
        print(X.shape)

        print("y shape before")
        print(y.shape)
        X = X[~na_mask,:]
        y = y[~na_mask]

        print("before shape")
        print(metadata_labels.shape)

        metadata_labels_temp = metadata_labels.loc[~na_mask,:]
        print(metadata_labels_temp.columns)

        if not bool_lodo:
            
            pc_scores = pc_table # get pc scores
            pc_scores = pc_scores[~na_mask,:]

            pc_scores_temp = pd.DataFrame(pc_scores,index = metadata_labels_temp.index)
            to_plot_correlation = pd.concat([metadata_labels_temp,pc_scores_temp],sort=False,ignore_index=True,axis=1)#,axis=0,ignore_index=True)
            to_plot_correlation.columns = np.append(np.array(metadata_labels_temp.columns),np.array(["PC" + str(i+1) for i in range(num_pcs)]))
            to_plot_correlation.to_csv(output_folders[d] + "pc_and_others.txt")
        
        
        print("Shape metadata temp")
        print(metadata_labels_temp.shape)
        # for each test train split in 5 fold cross validation
        
        train_it = 0


        # if doing lodo, designature groups for leave one dataset out
        if bool_lodo:
            print("lodo time")

            groups = np.array(metadata_labels[lodo_group])
            
            splitter = logo.split(X, y, groups)

        else:
            splitter = rskf.split(X, y)

        print("Compare number of splits")
        print("stratified")
        print(rskf.get_n_splits(X,y))
        print("stratified")
        groups = np.array(metadata_labels[lodo_group])
        print(logo.get_n_splits(X, y, groups))

        for train_index, test_index in splitter:

            # only run on the iteration in this grid cell (for paralle work)
            if train_it == train_it_input:
            
                test_train_start = timer()
                
                X_train, X_test = X[train_index,], X[test_index,]
                y_train, y_test = y[train_index], y[test_index]

                #print(metadata_labels_temp["sampleID"][test_index[0:50]])

                if bool_lodo:
                    pca_fit = pca.fit(X_train)
                    pc_scores_train = pca.transform(X_train)
                    pc_scores_test = pca.transform(X_test)

                else:
                    pc_scores_train, pc_scores_test =  pc_scores[train_index], pc_scores[test_index]

                if use_validation:

                    n_splits = 5
                    n_repeats = 1

                    X_train, X_val, y_train, y_val, pc_scores_train, pc_scores_val = train_test_split(X_train, y_train,pc_scores_train, test_size=0.30, random_state=1) 
          
                
                results_dict["number samples"] = []
                results_dict["number samples"].append(X_train.shape[0])
                # for each PC we regress out 
                for p in range(num_pcs): #range(3):#
                    results_dict["PC" + str(p)] = dict()
                       
                    if p == 0:
                        X_train_corrected = X_train
                        X_test_corrected = X_test
                        if use_validation:
                            X_val_corrected = X_val
                    else:
                        X_train_corrected = pca_regression(X_train,pc_scores_train[:,0:p])
                        X_test_corrected = pca_regression(X_test,pc_scores_test[:,0:p])
                        if use_validation:
                            X_val_corrected  = pca_regression(X_val,pc_scores_val[:,0:p])
                    
                    # perform grid search on train
                    best_train_model, best_params = RF_grid_search(X_train_corrected, y_train,parameter_dict)
                    
                    pickle.dump(best_train_model , open( output_folders[d] + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl", "wb" ) )
                     
                    
                test_train_end = timer()
                print("Finished one test train split")
                print(test_train_end - test_train_start)
            train_it += 1
    elif perform_MINERVA == 2: # predict with PCs directly
        # get PC scores
        pca = PCA(n_components=num_pcs,svd_solver='randomized')
        if not bool_lodo:
            # do the PCA thing
            temp = feature_table.transpose()
            pca.fit(temp)
            print("Time for pca cal")
            pca_start= timer()
            print("pca start clock")
            print(pca_start)
            pc_table = pca.transform(temp)  
            print(timer() - pca_start )
        feature_table_np = np.array(feature_table)
        labels_np = np.array(labels)
         
        dataset_start= timer()
        results_dict = dict()
        X = feature_table_np.transpose()
        y = labels_np
        
        na_mask = pd.isna(y)
        print("Xshape before")
        print(X.shape)

        print("y shape before")
        print(y.shape)
        X = X[~na_mask,:]
        y = y[~na_mask]

        print("before shape")
        print(metadata_labels.shape)

        metadata_labels_temp = metadata_labels.loc[~na_mask,:]
        print(metadata_labels_temp.columns)

        if not bool_lodo:
            
            pc_scores = pc_table # get pc scores
            pc_scores = pc_scores[~na_mask,:]

            pc_scores_temp = pd.DataFrame(pc_scores,index = metadata_labels_temp.index)
            to_plot_correlation = pd.concat([metadata_labels_temp,pc_scores_temp],sort=False,ignore_index=True,axis=1)#,axis=0,ignore_index=True)
            to_plot_correlation.columns = np.append(np.array(metadata_labels_temp.columns),np.array(["PC" + str(i+1) for i in range(num_pcs)]))
            to_plot_correlation.to_csv(output_folders[d] + "pc_and_others.txt")
        
        
        print("Shape metadata temp")
        print(metadata_labels_temp.shape)
        # for each test train split in 5 fold cross validation
        
        train_it = 0


        # if doing lodo, designature groups for leave one dataset out
        if bool_lodo:
            print("lodo time")

            groups = np.array(metadata_labels[lodo_group])
            
            splitter = logo.split(X, y, groups)

        else:
            splitter = rskf.split(X, y)

        print("Compare number of splits")
        print("stratified")
        print(rskf.get_n_splits(X,y))
        print("stratified")
        #groups = np.array(metadata_labels[lodo_group])
        #print(logo.get_n_splits(X, y, groups))

        for train_index, test_index in splitter:

            # only run on the iteration in this grid cell (for paralle work)
            if train_it == train_it_input:
            
                test_train_start = timer()
                
                X_train, X_test = X[train_index,], X[test_index,]
                y_train, y_test = y[train_index], y[test_index]

                #print(metadata_labels_temp["sampleID"][test_index[0:50]])

                if bool_lodo:
                    pca_fit = pca.fit(X_train)
                    pc_scores_train = pca.transform(X_train)
                    pc_scores_test = pca.transform(X_test)

                else:
                    pc_scores_train, pc_scores_test =  pc_scores[train_index], pc_scores[test_index]

                if use_validation:

                    n_splits = 5
                    n_repeats = 1

                    X_train, X_val, y_train, y_val, pc_scores_train, pc_scores_val = train_test_split(X_train, y_train,pc_scores_train, test_size=0.30, random_state=1) 
          
                
                results_dict["number samples"] = []
                results_dict["number samples"].append(X_train.shape[0])


                best_train_model, best_params = RF_grid_search(pc_scores_train, y_train,parameter_dict)
                    
                pickle.dump(best_train_model , open( output_folders[d] + special_name + file_output_string + "_PC" + str(num_pcs) + "_grid.pkl", "wb" ) )
                     
                    
                test_train_end = timer()
                print("Finished one test train split")
                print(test_train_end - test_train_start)

            train_it += 1
            
                

    all_datasets_dict["dataset" + str(d)] = results_dict
    
    dataset_end = timer()
    print("Finished one dataset")
    print(dataset_end - dataset_start)

# ...
end = timer()
print(end - start)

print("While time")
print(timer()- whole_time_start )   
#pickle.dump(all_datasets_dict , open( metadata_folder + special_name + file_output_string  + "_MINERVA_tt_grid.pkl", "wb" ) )
            
        
    


