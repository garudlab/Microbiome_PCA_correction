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
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.linear_model import ElasticNet,LinearRegression,ElasticNetCV
from sklearn import model_selection 
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
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

perform_MINERVA = int(args[11])

train_it_input = int(args[12])

perform_enet = bool(int(args[13]))
alpha_input = float(args[14])
l1ratio_input = float(args[15])

lodo_group = args[16]

if perform_enet:
    file_output_string  = "PREDenet_alpha" + str(alpha_input) +  "_l1ratio" + str(l1ratio_input) + "_trainit" + str(train_it_input)
else:
    file_output_string  = "PRED" + "_trainit" + str(train_it_input)


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

def pred_enet_grid_search(data,labels,param_dict):
    clf = ElasticNetCV(random_state=0,alphas = param_dict['alpha'],l1_ratio = param_dict['l1ratio'])
    clf.fit(data, labels)
    print("Best parameters set found on development set:")
    best_params = clf.get_params
    return clf,best_params


def pred_cv(data,labels):
    clf =  LinearRegression().fit(data, labels)
    return clf

def pred_enet_cv(data,labels,param_dict):
    clf = ElasticNet(random_state=0,alpha = param_dict['alpha_input'],l1_ratio = param_dict['l1ratio_input'])
    results = reg.score(X=data, y=labels)
    print(results)
    return(clf)


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
print(metadata[column_of_interest])

print(metadata[column_of_interest][0:5])
metadata[column_of_interest] = [float('Nan') if i == 'not applicable' or i == 'not provided' or 
                              i == "Not provided" or i== "Unspecified" or i == "Not applicable" 
                              else float(i) for i in list(metadata[column_of_interest])]

print(metadata[column_of_interest][0:5])
# temp_labels = np.array(metadata[column_of_interest])
# temp_labels = temp_labels.astype(np.float)

# metadata[column_of_interest] = temp_labels


non_nan_samples = metadata.index[np.invert(np.isnan(metadata[column_of_interest]))]


# Random forest stuff
n_splits = 5
n_repeats = n_repeats_input
import random
print ("Random number with seed 30")
random.seed(30)

rskf = model_selection.RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)

parameter_dict = {'alpha':[alpha_input],'l1ratio': [l1ratio_input]}


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

    if perform_MINERVA == 0 or perform_MINERVA == 3 or perform_MINERVA == 4:
        feature_table_np = np.array(feature_table)
        labels_np = np.array(labels)
        dataset_start= timer()
        results_dict = dict()
        X = feature_table_np.transpose()
        y = labels_np
        na_mask = pd.isna(y)
        
        X = X[~na_mask,:]
        y = y[~na_mask]
        metadata_labels_temp = metadata_labels.loc[~na_mask,:]
        groups = np.array(metadata_labels_temp[lodo_group])
        groups_one_hot = pd.get_dummies(groups,drop_first=True)

        ## do DCC
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

        # for each test train split in 5 fold cross validation
        train_it = 0
        for train_index, test_index in rskf.split(X, y):   
            #print("len train index")
            #print(train_index[0:5]) 
            #print(len(train_index))
            #print(len(test_index))

            if train_it == train_it_input:  #only run for desired training iteration
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

                print("not e et")
                print("Shape xtrain")
                print(X_train.shape)
                print("Shape ytrain")
                print(y_train.shape)


                if perform_enet:
                    best_train_model, best_params = pred_enet_grid_search(X_train, y_train,parameter_dict)


                else:
                    best_train_model = pred_cv(X_train, y_train)

                print("best_train_model_results")
                print("true labels")
                print(y_train[0:5])
                print("pred_results")
                pred_results = best_train_model.predict(X_train)
                print(pred_results[0:5])
                print( np.corrcoef(x=list(y_train),y=list(pred_results)))
                print("score")
                print(best_train_model.score(X_train,y_train))




                pickle.dump(best_train_model, open( output_folders[d] + special_name + file_output_string  + "_grid.pkl", "wb" ) )
                test_train_end = timer()
                print("Finished one test train split")
                print(test_train_end - test_train_start)
                
            train_it += 1
            



    elif perform_MINERVA == 1:
        # get PC scores
        pca = PCA(n_components=num_pcs,svd_solver='randomized')
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
        pc_scores = pc_table # get pc scores
        na_mask = pd.isna(y)
        
        X = X[~na_mask,:]
        y = y[~na_mask]
        pc_scores = pc_scores[~na_mask,:]
        # for each test train split in 5 fold cross validation
        
        train_it = 0
        for train_index, test_index in rskf.split(X, y):
            #print("train index")
            #print(train_index[0:5]) 

            if train_it == train_it_input:
            
                test_train_start = timer()
                #print(train_index)
                #print(test_index)
                
                X_train, X_test = X[train_index,], X[test_index,]
                y_train, y_test = y[train_index], y[test_index]
                pc_scores_train, pc_scores_test =  pc_scores[train_index], pc_scores[test_index]

                if use_validation:

                    n_splits = 5
                    n_repeats = 1

                    X_train, X_val, y_train, y_val, pc_scores_train, pc_scores_val = train_test_split(X_train, y_train,pc_scores_train, test_size=0.30, random_state=1) 
                    #print("First 5 xval")
                    #print(X_val[0:5])
                
                
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

                    if perform_enet:
                        best_train_model, best_params = pred_enet_grid_search(X_train_corrected, y_train,parameter_dict)
                    else:
                        best_train_model = pred_cv(X_train_corrected, y_train)
                    
                    pickle.dump(best_train_model , open( output_folders[d] + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl", "wb" ) )
                     
                    
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
            
        
    


