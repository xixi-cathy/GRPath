print("************************************* Detecting noteworthy cell types *************************************")
import numpy as np
import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import resample
from sklearn.metrics import roc_curve
import warnings
warnings.filterwarnings("ignore")
init_path = os.getcwd()+ '/../demo_data'
os.chdir(init_path)

data_name = "single_cell_data.txt"
label_name = "single_cell_label.txt"
iter_num = 3

def calculate_auc(X, y, cell_type, gene, sample_size = 200, sample_round = 10, fold = 4):
    X_temp = np.array(X.loc[gene, y["cell_type"]==cell_type])
    y_temp = np.array(y.loc[y["cell_type"]==cell_type, "label"])
    aucs = []
    for i in range(sample_round):
        # subsampling
        X_sample1, y_sample1 = resample(X_temp[y_temp==np.unique(y_temp)[1]], y_temp[y_temp==np.unique(y_temp)[1]], replace = True, n_samples = sample_size)
        X_sample0, y_sample0 = resample(X_temp[y_temp==np.unique(y_temp)[0]], y_temp[y_temp==np.unique(y_temp)[0]], replace = True, n_samples = sample_size)
        X_sample = np.append(X_sample0, X_sample1)
        X_sample = X_sample.reshape((len(X_sample), 1))
        y_sample = np.append(y_sample0, y_sample1)

        # 
        skf = StratifiedKFold(n_splits=fold)
        for train, test in skf.split(X_sample, y_sample):
            #clf = LogisticRegression(solver="liblinear",random_state=0)
            clf = RandomForestClassifier(max_depth=2, n_estimators = 10)
            clf.fit(X_sample[train], y_sample[train])
            fpr, tpr, threshholds = roc_curve(y_sample[test], clf.predict_proba(X_sample[test])[:, 1], pos_label = "cHF")
            a = metrics.auc(fpr, tpr)
            aucs.append(a)
    return np.mean(aucs)

X = pd.read_csv(data_name, sep = '\t')
y = pd.read_csv(label_name, sep = '\t')

alter_pathogenic = pd.read_csv("../alter_pathogenic_region_variant_egene.csv", sep=',')
alter_egenes = list(set(alter_pathogenic["heart_egenes"]))
alter_egenes.sort(key=list(alter_pathogenic["heart_egenes"]).index)
ref_pathogenic = pd.read_csv("../ref_pathogenic_region_variant_egene.csv", sep=',')
ref_egenes = list(set(ref_pathogenic["heart_egenes"]))
ref_egenes.sort(key=list(ref_pathogenic["heart_egenes"]).index)

egenes = alter_egenes + ref_egenes
genes = list(set(egenes).intersection(set(list(X.index))))
genes.sort(key=egenes.index)

celltypes = np.unique(y["cell_type"])
auc_summary = np.ndarray(shape=(len(genes), len(celltypes), iter_num))
for time in range(iter_num):
    print("Iteration: " + str(time+1))
    for i in range(len(celltypes)):
        print('\tCell type: '  + celltypes[i])
        for j in range(len(genes)):
            print('\t\tGene: '  + genes[j])
            auc_summary[j, i, time] = calculate_auc(X, y, celltypes[i], genes[j])
auc_mean = np.mean(auc_summary, axis=2)
auc_std = np.std(auc_summary, ddof=1, axis=2)

df_mean = pd.DataFrame(auc_mean)
df_mean.index = genes
df_mean.columns = celltypes
df_std = pd.DataFrame(auc_std)
df_std.index = genes
df_std.columns = celltypes
df_mean.to_csv("../processing/cell_type_auc_mean.csv")
df_std.to_csv("../processing/cell_type_auc_std.csv")
print("Done!")

