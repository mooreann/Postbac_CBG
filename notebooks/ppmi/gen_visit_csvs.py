import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

WRKDIR = '/labseq/projects/ppmi/wb_rna/plink/results/'
visits_oi = ['BL','V02','V04','V06','V08']

final_indexes = []

#opening index pair files, loading as list
for x in visits_oi:  
    final_index_file = ('/home/mooreank/index_beta_plots/index_pairs_{}.txt').format(x)
    with open(final_index_file) as f:
        finalindexlist = [line.rstrip('\n') for line in open(final_index_file)]
        for item in finalindexlist:
            final_indexes.append(item)

final_indexes = set(final_indexes)    


##run for each visit

test = 'BL'
#make dataframe of all data from each visit
filename_test = '{}eqtl.{}.fdr.plink.glm.linear'.format(WRKDIR,test)
test_all = pd.read_csv(filename_test, sep = '\s+')
test_all['INDEX PAIR'] = test_all['GENEID'].map(str)+"_"+ test_all['ID']
test_all['ZSCORE'] = test_all['BETA']/test_all['SE']
test_all['VISIT'] = test

test_sub_df = test_all[test_all['INDEX PAIR'].isin(final_indexes)]

test_filtered = pd.DataFrame()
for hybrid in (test_sub_df['INDEX PAIR']).tolist():
        #create sub df of just one hybrids rows
        hybrid_df = test_sub_df.loc[test_sub_df['INDEX PAIR'] == hybrid]
        #find min p value
        min_value = min(hybrid_df['P'])
        p_df = (hybrid_df.loc[hybrid_df['P'] == min_value])
        #find max absolute value beta
        p_df['ABS_BETA'] = abs(p_df['BETA'])
        p_df_sort = p_df.sort_values(by=['ABS_BETA'], ascending = False)
        index_beta = p_df_sort['BETA'].values[0]
        final_row = test_sub_df.loc[test_sub_df['BETA']==index_beta]
        test_filtered = test_filtered.append(final_row)
        
test_filtered.to_csv('/home/mooreank/BL_filtered_index.csv', index=False)









test = 'V02'
#make dataframe of all data from each visit
filename_test = '{}eqtl.{}.fdr.plink.glm.linear'.format(WRKDIR,test)
test_all = pd.read_csv(filename_test, sep = '\s+')
test_all['INDEX PAIR'] = test_all['GENEID'].map(str)+"_"+ test_all['ID']
test_all['ZSCORE'] = test_all['BETA']/test_all['SE']
test_all['VISIT'] = test
#dataframe of only index pairs

test_sub_df = test_all[test_all['INDEX PAIR'].isin(final_indexes)]

test_filtered = pd.DataFrame()
for hybrid in (test_sub_df['INDEX PAIR']).tolist():
        #create sub df of just one hybrids rows
        hybrid_df = test_sub_df.loc[test_sub_df['INDEX PAIR'] == hybrid]
        #find min p value
        min_value = min(hybrid_df['P'])
        p_df = (hybrid_df.loc[hybrid_df['P'] == min_value])
        #find max absolute value beta
        p_df['ABS_BETA'] = abs(p_df['BETA'])
        p_df_sort = p_df.sort_values(by=['ABS_BETA'], ascending = False)
        index_beta = p_df_sort['BETA'].values[0]
        final_row = test_sub_df.loc[test_sub_df['BETA']==index_beta]
        test_filtered = test_filtered.append(final_row)
        
test_filtered.to_csv('/home/mooreank/V02_filtered_index.csv', index=False)











test = 'V04'
#make dataframe of all data from each visit
filename_test = '{}eqtl.{}.fdr.plink.glm.linear'.format(WRKDIR,test)
test_all = pd.read_csv(filename_test, sep = '\s+')
test_all['INDEX PAIR'] = test_all['GENEID'].map(str)+"_"+ test_all['ID']
test_all['ZSCORE'] = test_all['BETA']/test_all['SE']
test_all['VISIT'] = test
#dataframe of only index pairs

test_sub_df = test_all[test_all['INDEX PAIR'].isin(final_indexes)]

test_filtered = pd.DataFrame()
for hybrid in (test_sub_df['INDEX PAIR']).tolist():
        #create sub df of just one hybrids rows
        hybrid_df = test_sub_df.loc[test_sub_df['INDEX PAIR'] == hybrid]
        #find min p value
        min_value = min(hybrid_df['P'])
        p_df = (hybrid_df.loc[hybrid_df['P'] == min_value])
        #find max absolute value beta
        p_df['ABS_BETA'] = abs(p_df['BETA'])
        p_df_sort = p_df.sort_values(by=['ABS_BETA'], ascending = False)
        index_beta = p_df_sort['BETA'].values[0]
        final_row = test_sub_df.loc[test_sub_df['BETA']==index_beta]
        test_filtered = test_filtered.append(final_row)
        
test_filtered.to_csv('/home/mooreank/V04_filtered_index.csv', index=False)















test = 'V06'
#make dataframe of all data from each visit
filename_test = '{}eqtl.{}.fdr.plink.glm.linear'.format(WRKDIR,test)
test_all = pd.read_csv(filename_test, sep = '\s+')
test_all['INDEX PAIR'] = test_all['GENEID'].map(str)+"_"+ test_all['ID']
test_all['ZSCORE'] = test_all['BETA']/test_all['SE']
test_all['VISIT'] = test
#dataframe of only index pairs

test_sub_df = test_all[test_all['INDEX PAIR'].isin(final_indexes)]

test_filtered = pd.DataFrame()
for hybrid in (test_sub_df['INDEX PAIR']).tolist():
        #create sub df of just one hybrids rows
        hybrid_df = test_sub_df.loc[test_sub_df['INDEX PAIR'] == hybrid]
        #find min p value
        min_value = min(hybrid_df['P'])
        p_df = (hybrid_df.loc[hybrid_df['P'] == min_value])
        #find max absolute value beta
        p_df['ABS_BETA'] = abs(p_df['BETA'])
        p_df_sort = p_df.sort_values(by=['ABS_BETA'], ascending = False)
        index_beta = p_df_sort['BETA'].values[0]
        final_row = test_sub_df.loc[test_sub_df['BETA']==index_beta]
        test_filtered = test_filtered.append(final_row)
        
test_filtered.to_csv('/home/mooreank/V06_filtered_index.csv', index=False)









test = 'V08'
#make dataframe of all data from each visit
filename_test = '{}eqtl.{}.fdr.plink.glm.linear'.format(WRKDIR,test)
test_all = pd.read_csv(filename_test, sep = '\s+')
test_all['INDEX PAIR'] = test_all['GENEID'].map(str)+"_"+ test_all['ID']
test_all['ZSCORE'] = test_all['BETA']/test_all['SE']
test_all['VISIT'] = test
#dataframe of only index pairs

test_sub_df = test_all[test_all['INDEX PAIR'].isin(final_indexes)]

test_filtered = pd.DataFrame()
for hybrid in (test_sub_df['INDEX PAIR']).tolist():
        #create sub df of just one hybrids rows
        hybrid_df = test_sub_df.loc[test_sub_df['INDEX PAIR'] == hybrid]
        #find min p value
        min_value = min(hybrid_df['P'])
        p_df = (hybrid_df.loc[hybrid_df['P'] == min_value])
        #find max absolute value beta
        p_df['ABS_BETA'] = abs(p_df['BETA'])
        p_df_sort = p_df.sort_values(by=['ABS_BETA'], ascending = False)
        index_beta = p_df_sort['BETA'].values[0]
        final_row = test_sub_df.loc[test_sub_df['BETA']==index_beta]
        test_filtered = test_filtered.append(final_row)
        
test_filtered.to_csv('/home/mooreank/V08_filtered_index.csv', index=False)
