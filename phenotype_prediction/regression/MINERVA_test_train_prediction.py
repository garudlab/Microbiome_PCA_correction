#!/usr/bin/env python
# coding: utf-8

# In[1]:


#args =  ["classifier.py","/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "CRC_k6&CRC_k7", "kmer", "BatchCorrected",  "bin_crc_normal",1,0,20,1,1] 

#python MINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc "CRC_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 20 1 1
# python MINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc "AGP_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 20 1 1

# # In[2]:
# from sklearn.ensemble import RandomForestClassifier
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
# already_trained_score = roc_auc_score(y_true = y_true, y_score = y_pred[:,1])

# already_trained_score
# reducing n_estimators reduced AUC on train, increasing min_samples_leaf reduced AUC on train
use_validation = True 



import sys
import pandas as pd
import utils
import numpy as np
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.ensemble import RandomForestClassifier
from sklearn import model_selection 
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from collections import Counter
from timeit import default_timer as timer
import pickle
import os 



args = sys.argv
print(args)


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

perform_enet = bool(int(args[12]))

if len(args) > 13:
    # really just means batch labels
    lodo_group = args[13]

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

def read_enet_grid_search_pc_version(folder,param_dict,train_it_input,pc,x,y):
    print("rweading")
    score_train = []
    models = dict()
    model_count = 0

    for alpha_input in param_dict['alpha']:
        for l1ratio_input in param_dict['l1ratio']:
            file_output_string  = "PREDenet_alpha" + str(alpha_input) +  "_l1ratio" + str(l1ratio_input) + "_trainit" + str(train_it_input)
            print(folder + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl")

            if os.path.isfile(folder + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl"):
                clf = pickle.load(open(folder + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl","rb"))
                already_trained_score = clf.score(x,y)
                print(already_trained_score)

                models[model_count] = clf
                score_train.append(already_trained_score)
                model_count += 1

                
            else:
                print("notfile")

    print(score_train)
    print(score_train.index(max(score_train)))
    best_model_index = score_train.index(max(score_train))
    best_model = models[best_model_index]
    best_params = best_model.get_params
    print("Best parameters set found on development set:")
    print(best_params)

    
    return best_model,best_params
def read_lin_pc_version(folder,param_dict,train_it_input,pc,x,y):
    print("rweading")

    file_output_string  = "PRED" + "_trainit" + str(train_it_input)
    print(folder + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl")

    if os.path.isfile(folder + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl"):
        clf = pickle.load(open(folder + special_name + file_output_string + "_PC" + str(p) + "_grid.pkl","rb"))
        already_trained_score = clf.score(x,y)
        print(already_trained_score)
        
    else:
        print("notfile")


    best_model = clf
    best_params = best_model.get_params
    print("Best parameters set found on development set:")
    print(best_params)

    
    return best_model,best_params

def read_lin(folder,train_it_input,x,y):



    file_output_string  = "PRED" + "_trainit" + str(train_it_input)
    print(folder + special_name + file_output_string + "_grid.pkl")
    if os.path.isfile(folder + special_name + file_output_string + "_grid.pkl"):
        clf = pickle.load(open(folder + special_name + file_output_string + "_grid.pkl","rb"))
        y_train_score = clf.score(x,y)
        already_trained_score = y_train_score
        
        print(already_trained_score)

    return clf

def read_enet_grid_search(folder,param_dict,train_it_input,x,y):
    print("rweading")
    score_train = []
    models = dict()
    model_count = 0
    for alpha_input in param_dict['alpha']:
        for l1ratio_input in param_dict['l1ratio']:
            file_output_string  = "PREDenet_alpha" + str(alpha_input) +  "_l1ratio" + str(l1ratio_input) + "_trainit" + str(train_it_input)
            print(folder + special_name + file_output_string + "_grid.pkl")

            if os.path.isfile(folder + special_name + file_output_string + "_grid.pkl"):
                clf = pickle.load(open(folder + special_name + file_output_string + "_grid.pkl","rb"))
                y_train_score = clf.score(x,y)
                already_trained_score = y_train_score
                print(already_trained_score)

                models[model_count] = clf
                score_train.append(already_trained_score)

                model_count += 1
            else:
                print("notfile")
def read_lin_model(folder,param_dict,train_it_input,x,y):
    print("rweading")

    file_output_string  = "PRED" + "_trainit" + str(train_it_input)
    print(folder + special_name + file_output_string + "_grid.pkl")

    if os.path.isfile(folder + special_name + file_output_string + "_grid.pkl"):
        clf = pickle.load(open(folder + special_name + file_output_string + "_grid.pkl","rb"))
        y_train_score = clf.score(x,y)
        already_trained_score = y_train_score
        print(already_trained_score)


    else:
        print("notfile")

                            



    best_model = clf
    best_params = clf.get_params
    print("Best parameters set found on development set:")
    print(best_params)

    
    return best_model,best_params

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


original_non_nan_samples = metadata.index[np.invert(np.isnan(metadata[column_of_interest]))]
# Random forest stuff
n_splits = 5
n_repeats = n_repeats_input
import random
print ("Random number with seed 30")
random.seed(30)

rskf = model_selection.RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)

parameter_dict = {'alpha':[ 0.025, 0.05, .125, .25, .5, 1., 2., 4.],'l1ratio':[0,.1, .5, .7, .9, .95, .99, 1]}


for d in range(len(study_names)): # range(1):#
    ### COPY INTO loop
    feature_table_dict = utils.load_feature_table([data_folders[d]],data_type = data_type)
    feature_table = feature_table_dict[0]

    print("feature table shape")
    print(feature_table.shape)

    

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
        


    #print(Counter(metadata[column_of_interest]))


    #########################################################################
    ###### COMMENTARY:  efining labels and binarize if not already     ######
    #########################################################################

    

    non_nan_samples = intersection(feature_table.columns,original_non_nan_samples)
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

        
        print("Shape original x")
        print(X.shape)
        print("Shape y")
        print(y.shape)

        metadata_labels_temp = metadata_labels.loc[~na_mask,:]
        groups = np.array(metadata_labels_temp[lodo_group])
        groups_one_hot = pd.get_dummies(groups,drop_first=True)

        ## perform DCC
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
            #print("train index")
            #print(train_index[0:5]) 

            if train_it == 0: 
                results_dict = dict()
                results_dict['train_best_params'] = dict()
                results_dict['train_score_trained'] = []
                
                results_dict['test_score_trained'] = []
                results_dict['val_score_trained']= []  

                results_dict['train_pearson_trained'] = []
                results_dict['val_pearson_trained'] = []  
                results_dict['test_pearson_trained'] = []  

            test_train_start = timer()
            X_train, X_test = X[train_index,], X[test_index,]
            y_train, y_test = y[train_index], y[test_index]
            if use_validation:
                n_splits = 5
                n_repeats = 1
                X_train, X_val, y_train, y_val = train_test_split(X_train, y_train,test_size=0.30, random_state=1) 
                #print("First 5 xval")
                #print(X_val[0:5])

            if train_it == 0: 
                results_dict["number samples"] = []
            results_dict["number samples"].append(X_train.shape[0])
            # for each PC we regress out 
            
                
            # perform grid search on train
            print("train Index" + str(train_it))
            if perform_enet:
                best_train_model, best_params = read_enet_grid_search(output_folders[d],parameter_dict, train_it, X_train, y_train)
            else:
                print("not e et")
                print("Shape xtrain")
                print(X_train.shape)
                print("Shape ytrain")
                print(y_train.shape)
                best_train_model, best_params = read_lin_model(output_folders[d],parameter_dict, train_it, X_train, y_train)


            print("best_train_model_results")
            print("true labels")
            print(y_train[0:5])
            print("pred_results")
            pred_results = best_train_model.predict(X_train)
            print(pred_results[0:5])
            print( np.corrcoef(x=list(y_train),y=list(pred_results)))
            print("score")
            print(best_train_model.score(X_train,y_train))

            # save best params
            results_dict['train_best_params'][train_it] = 1#best_params
            print("finished grid search: " + "train it " + str(train_it) )
            # get predictions on trained model
            
            already_trained_score = best_train_model.score(X_train, y_train)
            results_dict['train_score_trained'].append(already_trained_score)
            print("trained_model train Rf " + str(already_trained_score))  

            y_train_pred_prob = best_train_model.predict(X_train)
            already_trained_pearson= np.corrcoef(x=list(y_train),y=list(y_train_pred_prob))[0,1]
            results_dict['train_pearson_trained'].append(already_trained_pearson)
            print("trained pearson " + str(already_trained_pearson))  



            # # get predictions on newly trained model
            # newly_trained_auc = RF_cv(X_train,y_train,best_params)
            # results_dict['mean_train_cv_auc'].append(newly_trained_auc)
            # print("newly trained mean trained mean Rf " + str(np.mean(newly_trained_auc)))
            
            # validation metrics
            if use_validation:
                already_trained_val_score = best_train_model.score(X_val,y_val)
                results_dict['val_score_trained'].append(already_trained_val_score)
                print("trained_model validation Rf " + str(already_trained_val_score)) 
                y_val_pred_prob = best_train_model.predict(X_val)

                val_pearson= np.corrcoef(x=list(y_val),y=list(y_val_pred_prob))[0,1]
                results_dict['val_pearson_trained'].append(val_pearson)
                print("val pearson " + str(val_pearson))  


            # test metrics
            # get predictions on trained model
            already_trained_test_score = best_train_model.score(X_test,y_test)
            results_dict['test_score_trained'].append(already_trained_test_score)
            print("trained_model test RF" + str(already_trained_test_score))

            y_test_pred_prob = best_train_model.predict(X_test)

            test_pearson= np.corrcoef(x=list(y_test),y=list(y_test_pred_prob))[0,1]
            results_dict['test_pearson_trained'].append(test_pearson)
            print("test pearson " + str(test_pearson))  

            # # get predictions on newly trained model
            # test_RF = RF_cv(X_test,y_test,best_params)
            # results_dict['mean_test_cv_auc'].append(test_RF)
            # print("newly trained test mean Rf " + str(np.mean(test_RF)))
            
            all_datasets_dict["dataset" + str(d)] = results_dict  
            print("print all datasets dict")
            print(all_datasets_dict)
            #pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_tt_grid.pkl", "wb" ) )

            if use_validation:
                if not perform_enet:
                    pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_grid_VAL.pkl", "wb" ) )
                else:
                    pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_enet_grid_VAL.pkl", "wb" ) )
                            

            else:


                if not perform_enet:    
                    pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_grid.pkl", "wb" ) )
                else:
                    pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_enet_grid.pkl", "wb" ) )
                            

                
            train_it += 1
            test_train_end = timer()
            print("Finished one test train split")
            print(test_train_end - test_train_start)



    elif perform_MINERVA == 1:
        # get PC scores
        pca = PCA(n_components=num_pcs,svd_solver='randomized')
        # do the PCA thing
        temp = feature_table.transpose()
        pca.fit(temp)
        pc_table = pca.transform(temp)  
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
        
        starting_index = 0
        train_it = 0

        for train_index, test_index in rskf.split(X, y):
            test_train_start = timer()

            if train_it >= starting_index:
                #print("train index")
                #print(train_index[0:5])
                
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
                
                if train_it == starting_index: 
                    results_dict["number samples"] = []
                results_dict["number samples"].append(X_train.shape[0])
                # for each PC we regress out 
                for p in range(num_pcs): #range(3):#
                    if train_it == starting_index: 
                        results_dict["PC" + str(p)] = dict()
                        results_dict["PC" + str(p)]['train_best_params'] = dict()
                        results_dict["PC" + str(p)]['train_auc_trained'] = []
                        results_dict["PC" + str(p)]['test_auc_trained'] = []
                        results_dict["PC" + str(p)]['val_auc_trained']= []

                        results_dict["PC" + str(p)]['train_pearson_trained'] = []
                        results_dict["PC" + str(p)]['val_pearson_trained'] = []  
                        results_dict["PC" + str(p)]['test_pearson_trained'] = []  
                       
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
                        best_train_model, best_params = best_train_model, best_params = read_enet_grid_search_pc_version(output_folders[d],parameter_dict,train_it_input = train_it,pc=p,x=X_train_corrected,y = y_train)
                    else:
                        print("not e et")
                        print("Shape xtrain")
                        print(X_train.shape)
                        print("Shape ytrain")

                        print(y_train.shape)

                        #read_lin_pc_version(folder,param_dict,train_it_input,pc,x,y):
                        best_train_model, best_params = read_lin_pc_version(output_folders[d],parameter_dict, train_it, p,X_train, y_train)
                                        
                    # save best params
                    results_dict["PC" + str(p)]['train_best_params'][train_it] = 1
                    print("finished grid search: " + "train it " + str(train_it) + ", PC" + str(p))
                    # get predictions on trained model
                    #y_train_pred_prob = best_train_model.predict_proba(X_train_corrected)
                    already_trained_score = best_train_model.score(X_train_corrected,y_train)
                    results_dict["PC" + str(p)]['train_auc_trained'].append(already_trained_score)
                    print("trained_model train Rf " + str(already_trained_score))  


                    trained_pred_prob = best_train_model.predict(X_train_corrected)   
                    train_pearson= np.corrcoef(x=list(y_train),y=list(trained_pred_prob))[0,1]
                    results_dict["PC" + str(p)]['train_pearson_trained'].append(train_pearson)
                    print("train pearson " + str(train_pearson))            

                    # # get predictions on newly trained model
                    # newly_trained_auc = RF_cv(X_train_corrected,y_train,best_params)
                    # results_dict["PC" + str(p)]['mean_train_cv_auc'].append(newly_trained_auc)
                    # print("newly trained mean trained mean Rf " + str(np.mean(newly_trained_auc)))
                    
                    # validation metrics
                    if use_validation:
                        #y_val_pred_prob = best_train_model.predict_proba(X_val_corrected)
                        already_trained_val_score = best_train_model.score(X_val,y_val)
                        results_dict["PC" + str(p)]['val_auc_trained'].append(already_trained_val_score)
                        print("trained_model validation Rf " + str(already_trained_val_score)) 


                        val_pred_prob = best_train_model.predict(X_val_corrected)   
                        val_pearson= np.corrcoef(x=list(y_val),y=list(val_pred_prob))[0,1]
                        results_dict["PC" + str(p)]['val_pearson_trained'].append(val_pearson)
                        print("val pearson " + str(val_pearson))    

                    # test metrics
                    # get predictions on trained model
                    #y_test_pred_prob = best_train_model.predict_proba(X_test_corrected)
                    already_trained_test_score = best_train_model.score(X_test,y_test)
                    results_dict["PC" + str(p)]['test_auc_trained'].append(already_trained_test_score)
                    print("trained_model test RF" + str(already_trained_test_score))

                    test_pred_prob = best_train_model.predict(X_test_corrected)   
                    test_pearson= np.corrcoef(x=list(y_test),y=list(test_pred_prob))[0,1]
                    results_dict["PC" + str(p)]['test_pearson_trained'].append(test_pearson)
                    print("test pearson " + str(test_pearson))    

                    # # get predictions on newly trained model
                    # test_RF = RF_cv(X_test_corrected,y_test,best_params)
                    # results_dict["PC" + str(p)]['mean_test_cv_auc'].append(test_RF)
                    # print("newly trained test mean Rf " + str(np.mean(test_RF)))
                    
                    
                    all_datasets_dict["dataset" + str(d)] = results_dict
                    print("print all datasets dict")
                    print(all_datasets_dict)
                    #pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_tt_grid.pkl", "wb" ) )

                    if use_validation:
                        if not perform_enet:
                            pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_grid_VAL.pkl", "wb" ) )
                        else:
                            pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_enet_grid_VAL.pkl", "wb" ) )
                                    
     
                    else:


                        if not perform_enet:    
                            pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_grid.pkl", "wb" ) )
                        else:
                            pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_enet_grid.pkl", "wb" ) )
                                    

                     
                    
                    
            train_it += 1
            test_train_end = timer()
            print("Finished one test train split")
            print(test_train_end - test_train_start)

    all_datasets_dict["dataset" + str(d)] = results_dict
    
    dataset_end = timer()
    print("Finished one dataset")
    print(dataset_end - dataset_start)

# ...
end = timer()
print(end - start)

if use_validation:
    if not perform_enet:
        pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_grid_VAL.pkl", "wb" ) )
    else:
        pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_enet_grid_VAL.pkl", "wb" ) )
                
     
else:


    if not perform_enet:    
        pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_grid.pkl", "wb" ) )
    else:
        pickle.dump(all_datasets_dict , open( metadata_folder +"_" + special_name + "_MINERVA_prediction_enet_grid.pkl", "wb" ) )
                
            
    


