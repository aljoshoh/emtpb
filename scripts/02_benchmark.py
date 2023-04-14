#!/usr/bin/env python
# coding: utf-8

# In[4]:


import sys
import os
import os.path
from pathlib import Path
#os.chdir(Path(os.path.realpath('__file__')).parents[1])
#sys.path.append(Path(os.path.realpath('__file__')).parents[1])
current_dir = os.path.abspath(os.path.dirname('__file__'))
current_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
os.chdir(current_dir)
sys.path.append(current_dir)
print(current_dir)

os.chdir('/vol/emtpb/emtpb')
sys.path.append('/vol/emtpb/emtpb')


# In[5]:


import pandas as pd
import numpy as np
from scipy import stats
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import make_scorer
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
#import matplotlib.pyplot as plt
from sklearn.pipeline import make_pipeline
import joblib
import json
import sys
import math
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
import itertools

from econml.dml import CausalForestDML
from sklearn.linear_model import LassoCV

#run = "run6"
score = ["","_gsva","_tan_088","_secrier_065"]
model_type = ["eln","grf"]
response_type = ["","auc"]
cancer_type = range(28)

# create all possible combinations
combinations = list(itertools.product(score, model_type, response_type, cancer_type))

# create a DataFrame for each combination
dfs = []
for comb in combinations:
    df = pd.DataFrame(list(comb)).T
    df.columns = ["score", "model_type", "response_type", "cancer_type"]
    dfs.append(df)

# concatenate all DataFrames into a single DataFrame
df = pd.concat(dfs, ignore_index=True)

# initialize a counter
counter = 0
# iterate over the rows of the DataFrame
for i, row in df.iterrows():
    # increment the counter by 1 every 28th row
    if i % len(range(28)) == 0:
        counter += 1
    # add the counter value to a new column
    df.loc[i, 'run'] = str(counter+6)

# save DataFrame df
df.to_csv("metadata/paper/benchmark_paper_exp1.csv", index=True)

verbose = False # verbose=T does not work for CE
example = False # "1159-GDSC2" # False
cv = True


# In[6]:


# use one parameter dict
if len(sys.argv)!=2:
    argument = ['-',1]
else:
    print('cmd entry:', sys.argv)
    argument = sys.argv
which_run = int(argument[1]) - 1
print("emtpb: Parallel run "+str(which_run))
### costum
#which_run = 0
###

# assign values for run
print("emtpb: Updated manual parallel run "+str(which_run))
params = df.iloc[which_run]

run = 'run'+params['run']
score = params['score']
model_type = params['model_type']
response_type = params['response_type']
cancer_type = params['cancer_type']

EMTscores = pd.read_csv("metadata/EMTscores"+str(score)+".csv")
cancertypes = np.concatenate((['PANCAN'], EMTscores['TCGA Desc'].unique()))
preds_dir = "metadata/"+run+"/predictions/"+cancertypes[cancer_type]+"/"
model_dir = "metadata/"+run+"/models/"+cancertypes[cancer_type]+"/"
os.makedirs(preds_dir, exist_ok=True)
os.makedirs(model_dir, exist_ok=True)

mut = pd.read_csv("metadata/matrix_mut_"+cancertypes[cancer_type]+".csv")
resp_full = pd.read_csv("metadata/matrix_resp"+response_type+".csv")
cols = resp_full.columns[resp_full.columns.str.contains('GDSC')]


# In[ ]:


# Functions
def t_test(x, vadj=7.25): # 7.25 for 5x5 cv
    x = x[~np.isnan(x)]
    n = len(x)
    stat = np.mean(x) / np.sqrt(np.var(x)/n*vadj)
    return 2 * stats.t.sf(np.abs(stat), n-1)


def run_model_cv(X, y, outer_seed, inner_seed, model_dir, preds_dir, names, which_index, outer_folds = 5, inner_folds = 5, repeats = 5, baseline = False, save_models = False, model_type = "eln", check_for_trained = True, min_samples = 25):
    
    # predictions save path and check if its exists already
    preds_filename = os.path.join(preds_dir, f"predictions_{model_type}_{which_index}_{baseline}.csv")
    if check_for_trained:
        if Path(preds_filename).exists():
            print("emtpb: iteration "+str(index)+" already exists... skipping...")
            return
    
    # remove interesting variable for baseline
    if baseline:
        X = np.delete(X, 0, axis=1)
    
    # initialize list of all preds
    all_all_predictions = []
    all_all_truth = []
    all_all_folds = []
   
    # Initialize names
    names = pd.DataFrame(names)
    names = pd.concat([names]*repeats, ignore_index = True)
    
    # check size of X,y and save placeholder
    if mat.shape[0] < min_samples: # minimum 25 samples needed (5 test samples in outer cv and 4 in inner)
        print("Amount of columns is only "+str(mat.shape[0])+"... skipping")
        names.to_csv(preds_filename, index=False)
        return
    
    # initialize the outer cross-validation
    outer_cv = KFold(n_splits=outer_folds, shuffle=True, random_state=outer_seed)

    # initialize the inner cross-validation
    inner_cv = KFold(n_splits=inner_folds, shuffle=True, random_state=inner_seed)

    # repeat the outer cross-validation 5 times
    for repeat in range(repeats):
        
        # initialize the list to store the predictions
        all_predictions = []
        all_truth = []
        all_folds = []

        # set the seed for the outer cross-validation
        outer_cv.random_state = outer_seed + repeat

        # loop through the splits of the outer cross-validation
        for i, (train_index, test_index) in enumerate(outer_cv.split(X)):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            
            # initialize the model
            
            # elastic net
            if model_type == "eln":
                model = ElasticNetCV(l1_ratio=[0.01, .1, .5, .9, .95, 1],
                                     cv=inner_cv,
                                     random_state=inner_seed,
                                     tol=0.01)
                

            # random forest
            if model_type == "rf":
                model = RandomForestRegressor(random_state=inner_seed)
                
                
            # neural network
            if model_type == "mlp":
                model = MLPRegressor(max_iter=1000,
                                     random_state=inner_seed)
                # Define the hyperparameters to search over
                
                hidden_layer_sizes = [10, 50, 100, 150, 200]
                n_layers = [2, 3 , 4]
                layer_permutations = list(itertools.chain(*[itertools.permutations(hidden_layer_sizes, n) for n in n_layers]))
                
                param_grid = {
                    'hidden_layer_sizes': layer_permutations
                    #'activation': ['relu'] # 'identity', 'logistic', 'tanh', 
                    #'learning_rate_init': [0.001]
                }
                grid_search = RandomizedSearchCV(model, param_grid, cv=inner_folds, n_iter = 10, random_state = i+repeat*35)
                grid_search.fit(X, y)
                model = grid_search.best_estimator_
                
                # Get the best hyperparameters
                print(grid_search.best_params_)
                
                
            ############################    
            if model_type in ["eln","rf","mlp"]:
                
                # fit the model on the training data
                model.fit(X_train, y_train)
                
                # make predictions on the test set
                y_pred = model.predict(X_test)
                
                # save the predictions
                all_predictions.append(y_pred)
                all_truth.append(y_test)
                all_folds.append([i]*len(y_test))
                
            # causal modelling
            if model_type in ["grf"]:

                basemodel = ElasticNetCV(l1_ratio=[0.01, .1, .5, .9, .95, 1],
                                         cv=inner_cv,
                                         random_state=inner_seed,
                                         tol=0.01)
                
                # fit the model on the training data
                np.random.seed(123)
                model = CausalForestDML(discrete_treatment=False, random_state=123, model_y=basemodel, model_t=basemodel)

                if cv:
                    #tmp = np.array(y) # np.random.permutation(
                    X_red = np.delete(X_train,(0), axis=1) #
                    X_red_test = np.delete(X_test,(0), axis=1) #

                    # fit model
                    model.fit(Y = y_train, T = X_train[:,0], X = X_red) #model.fit(tmp, mat[:,0], X=matt)
                    summary = model.const_marginal_ate_inference(X_red_test) #, model.const_marginal_ate_interval(np.delete(mat,(0), axis=1))
                    #print(summary)
                    #print(summary.pvalue())
                    
                    # make predictions on the test set
                    #y_pred = model.predict(X_test)
                    y_pred = model.effect(X_red_test)

                    # save the predictions
                    all_predictions.append(y_pred)
                    all_truth.append(y_test)
                    all_folds.append([np.array2string(np.array([summary.pvalue(),model.const_marginal_ate(X_red_test),model.const_marginal_ate_interval(X_red_test)], dtype=object))]*len(y_test)) # transformback with: transformed_array = np.fromstring(my_array.replace('\n', '').replace(' ', ',').replace(',,',',').replace('(','').replace(')','')[1:-1], sep=',')
                    #all_folds.append([i]*len(y_test))
                    
                    pass
                    
                else:
                    # ignore train and test because we don't need it here
                    #tmp = np.array(y) # np.random.permutation(
                    X_red = np.delete(X,(0), axis=1) #

                    # fit model
                    model.fit(y, X[:,0], X=X_red) #model.fit(tmp, mat[:,0], X=matt)
                    summary = model.const_marginal_ate_inference(X_red) #, model.const_marginal_ate_interval(np.delete(mat,(0), axis=1))
                    #print(summary)
                    #print(summary.pvalue())
                    
                    # make predictions on the test set
                    #y_pred = model.predict(X_test)
                    y_pred = model.effect(X_red)

                    # save the predictions
                    all_predictions.append(y_pred)
                    all_truth.append(y)
                    all_folds.append([np.array2string(np.array([summary.pvalue(),model.const_marginal_ate(X_red),model.const_marginal_ate_interval(X_red)], dtype=object))]*len(y)) # transformback with: transformed_array = np.fromstring(my_array.replace('\n', '').replace(' ', ',').replace(',,',',').replace('(','').replace(')','')[1:-1], sep=',')
                    #all_folds.append([i]*len(y_test))
                    
                    break
                
                if cv:
                    pass
                else:
                    break # exit the loop after causal estimation without cv
                
            # save the predictions
            #all_predictions.append(y_pred)
            #all_truth.append(y_test)
            #all_folds.append([i]*len(y_test))

            # save the model
            model_filename = os.path.join(model_dir, f"model_{model_type}_{which_index}_{baseline}_{repeat}_{i}.joblib")
            if save_models:
                joblib.dump(model, model_filename)

        # concatenate all the predictions from cv
        all_predictions = np.concatenate(all_predictions)
        all_truth = np.concatenate(all_truth)
        all_folds = np.concatenate(all_folds)
        
        # save the cv results and append to resampling
        all_all_predictions.append(all_predictions)
        all_all_truth.append(all_truth)
        all_all_folds.append(all_folds)
    
    # concatenate all the predictions
    all_all_predictions = np.concatenate(all_all_predictions)
    all_all_truth = np.concatenate(all_all_truth)
    all_all_folds = np.concatenate(all_all_folds)
    
    names["preds"] = all_all_predictions
    names["truth"] = all_all_truth#np.concatenate([y]*repeats)
    names["repeat"] = np.concatenate([[i]*len(y) for i in range(repeats)])
    names["folds"] = all_all_folds
    names["repeatfold"] = names["repeat"].astype(str) + names["folds"].astype(str)
    
    # save predictions 
    names.to_csv(preds_filename, index=False)
    return names


# In[ ]:


# benchmark
for index in range(len(cols)):
    
    ### costum
    #index = 0
    ###
    if not not example:
        index = np.where(cols == "1559-GDSC2")[0][0]
    
    print("emtpb: modelling iteration "+str(index)+"...")
    #index = np.where(cols == "1559-GDSC2")[0][0]

    # get one response column
    resp = resp_full.copy()
    resp = resp[["COSMIC ID","TCGA Desc",cols[index]]]
    mat = pd.merge(EMTscores, resp, how='outer')
    if cancertypes[cancer_type] != "PANCAN":
        mat = mat[mat['TCGA Desc'] == cancertypes[cancer_type]]
    mat = pd.merge(mat, mut, how='outer')
    mat = pd.get_dummies(mat)
    mat_df = mat
    print("empb: shape before removing nan is "+str(mat.shape))
    mat.dropna(inplace=True)
    print("empb: shape after removing nan is "+str(mat.shape))
    if mat.shape[0] == 0:
        df = pd.DataFrame([[0]*mat.shape[1]],columns=mat.columns)
        mat = mat.append(df, ignore_index=True)
    names = mat['COSMIC ID']
    mat.drop(columns=['COSMIC ID'], inplace=True)
    response = mat.iloc[:,1] # mat.iloc[:,1] (drug response), # mat.iloc[:,0] (emt score)
    mat = mat.drop(mat.columns[1], axis=1)
    scale = StandardScaler()
    mat = scale.fit_transform(mat)

    for baseline in [True,False]: # [True, False]
        if model_type == "grf" and baseline == True:
            continue # skip baseline calcs for causal modelling
        df = run_model_cv(
            X = mat,
            y = np.array(response),
            outer_seed = 53+53, 
            inner_seed = 53,
            model_dir = model_dir,
            preds_dir = preds_dir,
            names = names,
            which_index = index,
            repeats = 5, #(1 if model_type == 'grf' else 5), # <------- SET TO 5
            baseline = baseline,
            model_type = model_type,
            save_models = False,
            check_for_trained = True,
            outer_folds = 5
        )

    if verbose:
        test = pd.read_csv(preds_dir+"predictions_"+str(model_type)+"_"+str(index)+"_False.csv")
        perfalse = np.array(test.groupby("repeatfold")[['preds','truth']].corr().iloc[1::2,0])
        perfalse[np.isnan(perfalse)] = 0
        print(perfalse)
        print(np.mean(perfalse)) # 
        test = pd.read_csv(preds_dir+"predictions_"+str(model_type)+"_"+str(index)+"_True.csv")
        pertrue = np.array(test.groupby("repeatfold")[['preds','truth']].corr().iloc[1::2,0])
        pertrue[np.isnan(pertrue)] = 0
        print(pertrue)
        print(np.mean(pertrue)) #
        print("t-test p="+str(t_test(perfalse-pertrue)))
        
    if not not example:
        break


# In[ ]:


print("emt_pb: done modelling "+cancertypes[cancer_types]+"for run"+str(which_run)+".")


# In[ ]:


#Check feature importances
#m = joblib.load("metadata/run1/models/SKCM/model_513_False_0_0.joblib")


# In[ ]:




