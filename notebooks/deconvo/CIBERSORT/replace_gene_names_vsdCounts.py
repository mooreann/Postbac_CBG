import pandas as pd
import numpy as np

WRKDIR = '/Users/mooreank/Desktop/Raph/Deconvolution/transcript_data/'

#load expression data
#quants_df = pd.read_hdf('ppmi.v121018.vst.quants.hdf5')
quants_df = pd.read_csv('{}vsdCountsPPMI_fullB38.tsv'.format(WRKDIR), sep='\t')


#load the gencode transcript annotations data
gencode = pd.read_csv('{}gencode.v29.transcripts.txt.gz'.format(WRKDIR), sep='\t', comment='#')
#gencode = pd.read_csv('gencode.v32.transcripts.txt', sep='\t', comment='#')


#create a dictionary of gene IDs (current feature/column labels) to gene names map
#id_dict = dict(zip(gencode['gene_id'], gencode['gene_name'])) 
id_dict = dict(zip(gencode['gene_id'], gencode['gene_name'])) 
#rename the features/columns using the dict (mapping)
quants_df.replace(id_dict, inplace=True)

quants_df.to_csv('genes.vsdCountsPPMI_fullB38.tsv', sep='\t')
