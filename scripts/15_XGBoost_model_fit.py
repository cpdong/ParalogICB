#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 8 2023
@author: cpdong
"""

import argparse, os;
import pandas as pd
import numpy as np
import pickle, time
from pathlib import Path
from itertools import product
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier
import shap


if __name__ == '__main__':
    cpus = 8
    familySize_cutoffs = 20;
    train_prop = 0.8; 
    n_tree = 400;
    seed = 149
    params = {'objective': 'reg:logistic', 'max_depth': 8, 'min_child_weight': 1, 'min_split_loss': 0,
              'subsample': 1, 'colsample_bytree': 0.8, 'learning_rate': 0.2, 'n_estimators': n_tree}

    dataAll = pd.read_csv('Model_Inputs.tsv', header=0, index_col=0, sep='\t')
    
    data_train_n_test = dataAll[(data_train_n_test['familysize']<= familySize_cutoffs) &
                                          (data_train_n_test['BLASTP_identity'] >= 30) &
                                          #(data_train_n_test['BLASTP_evalue'] >= 5) & 
                                          (data_train_n_test['Rmsd'].notnull()) &
                                          (data_train_n_test['TCGA_surv_SKCM'].notnull()) &
                                          (data_train_n_test['DepMap_synthetic_All_pval'].notnull()) &
                                          (data_train_n_test['paraTKO'] >=0)]

    data_train_n_test = data_train_n_test.fillna(0)
    #data_train_n_test.to_csv('data_train_n_test_data.tsv', index=True, sep='\t')
    features = data_train_n_test.columns[:-1]
    
    train, test = df2[df2['is_train']==True], df2[df2['is_train']==False]
    y_train = pd.factorize(train['paraTKO'], sort=True)[0]
    x_train = train[features]

    # internal validation - 1
    y_test = pd.factorize(test['paraTKO'], sort=True)[0]
    x_test = test[features]

    estimator = XGBClassifier(n_jobs = cpus, random_state = 8).set_params(**params)    
    estimator.fit(x_train, y_train); # fittinga
    feat_imp = pd.DataFrame({'feature_names': estimator.feature_names_in_,
                                'feature_importances': estimator.feature_importances_})
    feat_imp.to_csv('model_XGB_feature_importance.tsv', index=False, sep='\t')

    x_test_pred = estimator.predict_proba(x_test)
    internal_auc = roc_auc_score(y_test, x_test_pred[:, 1])
    print('internal AUC:', internal_auc)
    test['proba'] = x_test_pred[:, 1]
    test.to_csv('internal_test_data_praba_result.tsv', index=False, sep='\t')
