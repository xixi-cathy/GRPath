import pandas as pd
import os
import sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy import stats
import math
import random
from multiprocessing import Pool
from functools import partial
from collections import Counter
import warnings
warnings.filterwarnings("ignore")

init_path = os.getcwd()+ '/../demo_data'
os.chdir(init_path)
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
snps_name = "GWAS.txt"
variant_openness_name = "variant_openness.txt"
variant_allele_name = "variant_allele.txt"
variant_pos_name = "variant_pos.txt"
RE_openness_name = "RE_openness.txt"
phenotype_name = "phenotype.txt"
chunksize = 1000000
kernels = 10
iters = 5

print("************************************* Processing matrices per region *************************************")
snps = pd.read_csv(snps_name, sep='\t')
rspos = pd.read_csv(variant_pos_name, sep='\t')
df_ref = pd.read_csv(RE_openness_name, sep='\t', header=None)
df_RE = df_ref.iloc[:, 0:2]
for j in range(snps.shape[0]):
    rs = snps.iloc[j, 0]
    chrome = 'chr' + str(snps.iloc[j, 1])
    position = snps.iloc[j, 2]
    
    df_score = pd.read_csv(variant_openness_name, sep='\t', chunksize=chunksize, header=None)
    result = pd.DataFrame()
    allele_result = pd.DataFrame()
    count = 0
    for chunk in df_score:
        allele = pd.read_csv(variant_allele_name, sep='\t', header=None, skiprows = count * chunksize, nrows=chunksize)
        allele.index = allele.index + count * chunksize
        chunk_new = pd.DataFrame()
        allele_new = pd.DataFrame()
        # 200kb region surrounding the GWAS SNP
        ind = chunk.loc[(chunk.iloc[:, 0] == chrome) & (chunk.iloc[:, 4]>=position-100000) & (chunk.iloc[:, 4]<=position+100000), :].index.tolist()
        chunk_new = chunk.loc[(chunk.iloc[:,  0] == chrome) & (chunk.iloc[:,  4]>=position-100000) & (chunk.iloc[:,  4]<=position+100000), :]
        allele_new = allele.loc[ind, :]
        result = result.append(chunk_new)
        allele_result = allele_result.append(allele_new)
        count = count + 1
    if result.shape[0]==0:
        print("No variant found within the 200kb region around " + rs + "!")
        continue
    result = result.reset_index(drop=True)
    allele_result = allele_result.reset_index(drop=True)
       
    var_RE = result.loc[:, 0:2].drop_duplicates()
    ref_result = df_ref.iloc[df_RE.reset_index().merge(var_RE, how="right").set_index('index').index.tolist(), :]
    var_RE.index = range(len(var_RE))
    ref_result.index = range(len(var_RE))
    if not os.path.exists(init_path + '/../processing'):
        os.mkdir(init_path + '/../processing')
    savepath = init_path + '/../processing'
    ref_result.to_csv(savepath + '/RE_openness_'+ rs + ".txt", index=False, sep='\t', header=False)
    
    df_allele_val = allele_result.loc[:, 8:]
    alter_allele = df_allele_val.apply(lambda x: x.sum(), axis=1)
    ind = result.loc[(alter_allele > 10) & (alter_allele < df_allele_val.shape[1]-10), :].index.tolist()
    result = result.loc[(alter_allele > 10) & (alter_allele < df_allele_val.shape[1]-10), :]
    allele_result = allele_result.loc[(alter_allele > 10) & (alter_allele < df_allele_val.shape[1]-10), :]
    if result.shape[0]==0:
        print("No variant left within the 200kb region around " + rs + "after filtering!")
        continue    
    result.to_csv(savepath + '/variant_openness_'+ rs + ".txt", index=False, sep='\t', header=False)
    allele_result.to_csv(savepath + '/variant_allele_'+ rs + ".txt", index=False, sep='\t', header=False)
    print("Number " + str(j+1) + ": " + rs + " done.")
print("Done!")

print("************************************* Deciding causal regions and variants *************************************")
os.chdir(init_path)
samples_all = pd.read_csv(phenotype_name, sep = '\t')
gwas_snps = list(snps.iloc[:,0])
def cal_pval():
    # Split data: 2/5 for beta, 3/5 for lambda and VCS score
    os.chdir(init_path + "/../processing")
    negative_samples = samples_all.loc[samples_all['Label']==0, :]
    negative_donors = negative_samples['Donor'].drop_duplicates()
    positive_samples = samples_all.loc[samples_all['Label']==1, :]
    positive_donors = positive_samples['Donor'].drop_duplicates()
    beta_neg_num = round(0.4*len(negative_donors))
    beta_pos_num = round(0.4*len(positive_donors))
    beta_neg_donors = pd.DataFrame(negative_donors[random.sample(list(negative_donors.index), beta_neg_num)])
    beta_neg_donors.columns = ['Donor']
    beta_pos_donors = pd.DataFrame(positive_donors[random.sample(list(positive_donors.index), beta_pos_num)])
    beta_pos_donors.columns = ['Donor']
    beta_samples_origin = pd.merge(left=samples_all, right=beta_neg_donors, on='Donor')
    beta_samples_origin = beta_samples_origin.append(pd.merge(left=samples_all, right=beta_pos_donors, on='Donor'))
        
    pool = Pool(kernels)
    partial_work = partial(cal_pval_per_snp, beta_samples_origin = beta_samples_origin)
    snps, min_ps_greater, p_variant_nums_greater, causal_variantss_greater, min_ps_less, p_variant_nums_less, causal_variantss_less = zip(*pool.map(partial_work, [snp for snp in gwas_snps]))
    pool.close()
    pool.join()
    results_greater = pd.DataFrame({'snp':snps, 'variant_num':p_variant_nums_greater, 'min_p':min_ps_greater, 'causal_variants':causal_variantss_greater})
    results_less = pd.DataFrame({'snp':snps, 'variant_num':p_variant_nums_less, 'min_p':min_ps_less, 'causal_variants':causal_variantss_less})

    return results_greater, results_less

def cal_pval_per_snp(snp, beta_samples_origin):
    print('\tCalculating ' + snp + " surrounding region...")
    RE_openness_name = "RE_openness_" + snp + ".txt"
    variant_openness_name = 'variant_openness_'+snp+'.txt'
    variant_allele_name = 'variant_allele_'+snp+'.txt'
    
    # Keep samples with label
    ind = samples_all['Sample_var'].drop_duplicates().tolist()
    ind = [i + 8 for i in ind]
    variant_openness = pd.read_csv(variant_openness_name, sep = '\t', header = None)
    variant_allele = pd.read_csv(variant_allele_name, sep = '\t', header = None)
    variant_openness = pd.concat([variant_openness.loc[:, 0:7], variant_openness.loc[:, ind]], axis=1)
    variant_allele = pd.concat([variant_allele.loc[:, 0:7], variant_allele.loc[:, ind]], axis=1)
    RE_openness = pd.read_csv(RE_openness_name, sep = '\t', header = None)
    ind = samples_all['Sample_RE'].tolist()
    ind = [i + 3 for i in ind]
    RE_openness = pd.concat([RE_openness.loc[:, 0:2], RE_openness.loc[:, ind]], axis=1)

    # Keep variants with #alternate allele>10
    alter_allele = variant_allele.loc[:, 8:].apply(lambda x: x.sum(), axis=1)
    variant_openness = variant_openness.loc[alter_allele>10, :]
    variant_allele = variant_allele.loc[alter_allele>10, :]
    var_RE = variant_openness.loc[:, 0:2].drop_duplicates()
    RE_openness_filtered = pd.DataFrame()
    for i in range(len(var_RE)):
        RE_openness_temp = RE_openness.loc[(RE_openness.iloc[:, 0] == var_RE.iloc[i, 0]) & (RE_openness.iloc[:, 1] == var_RE.iloc[i, 1]) 
                               & (RE_openness.iloc[:, 2] == var_RE.iloc[i, 2]), :]
        RE_openness_filtered = RE_openness_filtered.append(RE_openness_temp, ignore_index=True)
    
    # Ready the RE_openenss matrix to regress beta
    cols = np.array(beta_samples_origin['Sample_RE'])+3
    beta_RE_openness = RE_openness_filtered.loc[:, cols]
    beta_RE_openness = pd.concat([RE_openness_filtered.loc[:, 0:2], beta_RE_openness], axis = 1)
    beta_label = pd.DataFrame(beta_samples_origin['Label'])

    # Ready the variant_openness matrix to calculate VCS score
    beta_samples = beta_samples_origin.drop_duplicates(['Donor'])
    samples = samples_all.drop_duplicates(['Donor'])
    lambda_samples = beta_samples.append(samples)
    lambda_samples = lambda_samples.drop_duplicates(['Donor'],keep=False)
    cols = np.array(lambda_samples['Sample_var'])+8
    lambda_variant_openness = variant_openness.loc[:, cols]
    lambda_variant_openness = pd.concat([variant_openness.loc[:, 0:7], lambda_variant_openness], axis = 1)
    lambda_variant_allele = variant_allele.loc[:, cols]
    lambda_variant_allele = pd.concat([variant_allele.loc[:, 0:7], lambda_variant_allele], axis = 1)
    lambda_label = pd.DataFrame(lambda_samples['Label'])

    df_X = beta_RE_openness
    X = np.array(df_X.iloc[:, 3:]).T
    df_label = beta_label
    y = np.array(df_label.iloc[:, 0])
    df_score = lambda_variant_openness
    df_allele = lambda_variant_allele
    
    # beta
    scaler = StandardScaler().fit(X)
    X_transformed = scaler.transform(X)
    beta = pd.DataFrame(LogisticRegression(C=0.1).fit(X_transformed, y).coef_.T)
    beta = pd.concat([df_X.loc[:, 0:2], beta], axis = 1)

    # lambda
    ## first, filter variants
    df_allele_filtered = df_allele.loc[np.var(df_score.iloc[:, 8:], axis=1)>1e-10, :]
    df_allele_filtered = df_allele_filtered.iloc[:, 8:]
    df_allele_filtered = df_allele_filtered.reset_index(drop=True)
    df_score_filtered = df_score.loc[np.var(df_score.iloc[:, 8:], axis=1)>1e-10, :]
    df_score_filtered = df_score_filtered.reset_index(drop=True)
    variant_info = df_score_filtered.loc[:, 0:7]
    df_score_filtered = df_score_filtered.iloc[:, 8:]
    lambda_ = []
    for i in range(df_score_filtered.shape[0]):
        temp_ind = np.array(df_allele_filtered)[i,:]==0
        ref_mean = np.mean(np.array(df_score_filtered)[i,temp_ind])
        temp_ind = np.array(df_allele_filtered)[i,:]==1
        alter_mean = np.mean(np.array(df_score_filtered)[i,temp_ind])
        lambda_i = abs(ref_mean-alter_mean)
        lambda_.append(lambda_i)
    lambda_ = pd.DataFrame(lambda_)
    lambda_ = pd.concat([variant_info, lambda_], axis = 1)

    # VCS score
    vcs = []
    for i in range(len(lambda_)):
        beta_i = beta.loc[(beta.iloc[:, 0]==lambda_.iloc[i, 0]) & (beta.iloc[:, 1]==lambda_.iloc[i, 1]) & (beta.iloc[:, 2]==lambda_.iloc[i, 2]), :]
        beta_i = beta_i.iloc[:, 3].iloc[0]
        vcs_i = abs(beta_i * lambda_.iloc[i, 8])
        vcs.append(vcs_i)
    vcs = pd.DataFrame(vcs)
    vcs = pd.concat([variant_info, vcs], axis = 1)

    # Rank the variants
    vcs.columns = [0, 1, 2, 3, 4, 5, 6, 7, 'vcs']
    vcs_descend = vcs.sort_values(by='vcs', ascending=False)
    df_allele_filtered_descend = df_allele_filtered.iloc[vcs_descend.index.tolist(), :]
    y = np.array(lambda_label.iloc[:, 0])
    vcs_descend["pos"] = vcs_descend[3].map(str) + "_" + vcs_descend[5].map(str) + "_" + vcs_descend[6].map(str) + "_"+ vcs_descend[7].map(str) + "_b37"
    vcs_descend = pd.merge(vcs_descend, rspos, how="left", on="pos")


    # Odds ratio
    num = int(vcs_descend.shape[0]/2)
    ## Select top variants and calculate the OR
    OR_top_list = []
    for rownum in range(num):
        ref_all = sum(df_allele_filtered_descend.iloc[rownum, ]==0)
        ref_disease = sum((df_allele_filtered_descend.iloc[rownum, ]==0) & (y==1))
        ref_normal = sum((df_allele_filtered_descend.iloc[rownum, ]==0) & (y==0))
        alter_all = sum(df_allele_filtered_descend.iloc[rownum, ]==1)
        alter_disease = sum((df_allele_filtered_descend.iloc[rownum, ]==1) & (y==1))
        alter_normal = sum((df_allele_filtered_descend.iloc[rownum, ]==1) & (y==0))
        stat = pd.DataFrame({'Disease':{'Reference allele': ref_disease,'Alternate allele':alter_disease},
                           'Normal':{'Reference allele': ref_normal,'Alternate allele':alter_normal}})
        if (ref_disease !=0) & (alter_normal != 0) & (ref_normal != 0):
            OR = (alter_disease/alter_normal)/(ref_disease/ref_normal)
            OR_top_list.append(OR)
    ## Select bottom variants and calculate the OR
    bottom_ind = range(vcs_descend.shape[0]-1, vcs_descend.shape[0]-1-num, -1)
    OR_bottom_list = []
    for rownum in bottom_ind:
        ref_all = sum(df_allele_filtered_descend.iloc[rownum, ]==0)
        ref_disease = sum((df_allele_filtered_descend.iloc[rownum, ]==0) & (y==1))
        ref_normal = sum((df_allele_filtered_descend.iloc[rownum, ]==0) & (y==0))
        alter_all = sum(df_allele_filtered_descend.iloc[rownum, ]==1)
        alter_disease = sum((df_allele_filtered_descend.iloc[rownum, ]==1) & (y==1))
        alter_normal = sum((df_allele_filtered_descend.iloc[rownum, ]==1) & (y==0))
        stat = pd.DataFrame({'Disease':{'Reference allele': ref_disease,'Alter allele':alter_disease},
                           'Normal':{'Reference allele': ref_normal,'Alter allele':alter_normal}})
        if (ref_disease !=0) & (alter_normal != 0) & (ref_normal != 0):
            OR = (alter_disease/alter_normal)/(ref_disease/ref_normal)
            OR_bottom_list.append(OR)

    variant_num = range(math.ceil(0.02*vcs_descend.shape[0]), num)  # Use 2%-50% variants
    causal_variants_greater = []
    min_p_greater = 1
    p_variant_num_greater = 0
    causal_variants_less = []
    min_p_less = 1
    p_variant_num_less = 0
    for j in variant_num:
        try:
            data = pd.DataFrame({'Top':OR_top_list[:j], 'Bottom':OR_bottom_list[:j]})
            ## Calculate p-value
            t_greater, p_greater = stats.mannwhitneyu(OR_top_list[:j], OR_bottom_list[:j], 
                                                      alternative="greater")
            t_less, p_less = stats.mannwhitneyu(OR_top_list[:j], OR_bottom_list[:j], 
                                                      alternative="less")
            if p_greater <= min_p_greater:
                min_p_greater = p_greater
                p_variant_num_greater = j
                #if p_greater<0.05:
                causal_variants_greater = list(set(vcs_descend["mapped_snp"][:j].dropna().tolist()))
            if p_less <= min_p_less:
                min_p_less = p_less
                p_variant_num_less = j
                #if p_less<0.05:
                causal_variants_less = list(set(vcs_descend["mapped_snp"][:j].dropna().tolist()))
        except:
            continue
    return snp, min_p_greater, p_variant_num_greater, causal_variants_greater, min_p_less, p_variant_num_less, causal_variants_less

greater = pd.DataFrame()
less = pd.DataFrame()
for t in range(iters):
    print("Iteration: " + str(t+1))
    temp_greater, temp_less = cal_pval()
    greater = greater.append(temp_greater)
    less = less.append(temp_less)

nums_greater = []
ps_greater = []
fisher_ps_greater = []
ccausal_variants_greater = []
nums_less = []
ps_less = []
fisher_ps_less = []
ccausal_variants_less = []
for snp in gwas_snps:
    df_greater_temp = greater.loc[greater["snp"]==snp,:]
    min_ps_greater = list(df_greater_temp["min_p"])
    variant_nums_greater = list(df_greater_temp["variant_num"])
    causal_variants_greater = list(df_greater_temp["causal_variants"])    
    greater_max_index = min_ps_greater.index(max(min_ps_greater))
    min_ps_greater.remove(max(min_ps_greater))
    variant_nums_greater.remove(variant_nums_greater[greater_max_index])
    causal_variants_greater.remove(causal_variants_greater[greater_max_index])
    greater_min_index = min_ps_greater.index(min(min_ps_greater))
    min_ps_greater.remove(min(min_ps_greater))
    variant_nums_greater.remove(variant_nums_greater[greater_min_index])
    causal_variants_greater.remove(causal_variants_greater[greater_min_index])
    min_p_greater = np.array(min_ps_greater).mean()
    ps_greater.append(min_p_greater)
    _, fisher_p_greater = stats.combine_pvalues(min_ps_greater, method='fisher')
    fisher_ps_greater.append(fisher_p_greater)
    variant_num_greater = np.array(variant_nums_greater).mean()
    nums_greater.append(variant_num_greater)
    variants_unlist = dict(Counter([v for t in causal_variants_greater for v in t]))
    causal_variants_greater = [key for key,value in variants_unlist.items() if value > 1]
    ccausal_variants_greater.append(causal_variants_greater)
    
    df_less_temp = less.loc[less["snp"]==snp,:]
    min_ps_less = list(df_less_temp["min_p"])
    variant_nums_less = list(df_less_temp["variant_num"])
    causal_variants_less = list(df_less_temp["causal_variants"])
    less_max_index = min_ps_less.index(max(min_ps_less))
    min_ps_less.remove(max(min_ps_less))
    variant_nums_less.remove(variant_nums_less[less_max_index])
    causal_variants_less.remove(causal_variants_less[less_max_index])
    less_min_index = min_ps_less.index(min(min_ps_less))
    min_ps_less.remove(min(min_ps_less))
    variant_nums_less.remove(variant_nums_less[less_min_index])
    causal_variants_less.remove(causal_variants_less[less_min_index])
    min_p_less = np.array(min_ps_less).mean()
    ps_less.append(min_p_less)
    _, fisher_p_less = stats.combine_pvalues(min_ps_less, method='fisher')
    fisher_ps_less.append(fisher_p_less)
    variant_num_less = np.array(variant_nums_less).mean()
    nums_less.append(variant_num_less)
    variants_unlist = dict(Counter([v for t in causal_variants_less for v in t]))
    causal_variants_less = [key for key,value in variants_unlist.items() if value > 1]
    ccausal_variants_less.append(causal_variants_less)
    
cccausal_variants_greater = []
cccausal_variants_less = []
for i in range(len(ccausal_variants_greater)):
    cccausal_variants_greater.append(" ".join(ccausal_variants_greater[i]))
    cccausal_variants_less.append(" ".join(ccausal_variants_less[i]))

output_greater = pd.DataFrame({'snp':gwas_snps, 'p_fisher':fisher_ps_greater, 'causal_variants':cccausal_variants_greater})
output_greater = output_greater.sort_values(by = ['p_fisher'], ascending=True)
output_less = pd.DataFrame({'snp':gwas_snps, 'p_fisher':fisher_ps_less, 'causal_variants':cccausal_variants_less})
output_less = output_less.sort_values(by = ['p_fisher'], ascending=True)
output_greater.to_csv(savepath + "/region_alter_pathogenic.csv", index=False, sep=',', header=True)
output_less.to_csv(savepath + "/region_reference_pathogenic.csv", index=False, sep=',', header=True)
print("Done!")

