visits_oi = ["BL","V02", "V04","V06","V08"]
WRKDIR = '/labseq/projects/ppmi/wb_rna/plink/results/'

import sys
visit = sys.argv[1]
#print(visit)

##upload BL file
import numpy as np
import pandas as pd

#filename = '/labseq/projects/ppmi/wb_rna/plink/results/eqtl.BL.oi.plink.glm.linear'
filename = '{}eqtl.{}.oi.plink.glm.linear'.format(WRKDIR,visit)
print(filename)

visit_df = pd.read_csv(filename, sep = '\s+')
    
##filter out pairs with BH over 1%
#visit_df = visit_df_unfilt.loc[visit_df_unfilt['bh_fdr'] <= 0.01]
##adding column of GENEID_ID
visit_df['hybrid'] = visit_df['GENEID'].map(str)+"_"+ visit_df['ID']
#visit_df['ZSCORE'] = visit_df['BETA']/visit_df['SE']

#
#get test set
gene_set = list(set(visit_df["GENEID"]))
#test_gene_set = gene_set[0:100]


index_pairs = []
#index_zscores = []
#index_fdrs = []
#
for gene in gene_set:
    #create sub df of just one gene's variants
    gene_df = visit_df.loc[visit_df['GENEID'] == gene]
    
    #find min p value
    min_value = min(gene_df['P'])
    p_df = (gene_df.loc[gene_df['P'] == min_value])

    #find max beta value
    p_df['ABS_BETA'] = abs(p_df['BETA'])
    p_df_sort = p_df.sort_values(by=['ABS_BETA'], ascending = False)

    index_pair= p_df_sort['hybrid'].values[0]
    #index_zscore = p_df_sort['ZSCORE'].values[0]
    #index_fdr = p_df_sort['bh_fdr'].values[0]
  
    index_pairs.append(index_pair)
    #index_zscores.append(index_zscore)
    #index_fdrs.append(index_fdr)


with open ("index_pairs_"+str(visit)+".txt", 'w') as f:
    for item in index_pairs:
        f.write("%s\n" % item)
        
#with open (str(visit)+"_index_zscores.txt", 'w') as f:
#    for item in index_zscores:
#        f.write("%s\n" % item)
#        
#with open (str(visit)+"_index_fdrs.txt", 'w') as f:
#    for item in index_fdrs:
#        f.write("%s\n" % item)
