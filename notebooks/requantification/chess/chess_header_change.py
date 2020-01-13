import pandas as pd
import numpy as np

WRKDIR = '/Users/mooreank/Desktop/Raph/requant/chess/'

chess = pd.read_csv('{}chess2.2_and_refseq.gtf'.format(WRKDIR), sep='\t', header = None)
chess = pd.DataFrame(chess)
chess.columns = ['Chrom','Source','Type','Start','End','','','','Split']

chess['Split'] = chess['Split'].str.rstrip(';')
chess[['Transcript_id', 'Gene_id','Gene_name']] = chess.Split.str.split(";",expand=True) 
del chess['Split']

chess.loc[chess['Gene_name'].isnull(), 'Gene_name'] = 'None'
chess['Transcript_id'] = chess['Transcript_id'].map(lambda x: x.lstrip('transcript_id "').rstrip('"'))
chess['Gene_id'] = chess['Gene_id'].map(lambda x: x.lstrip('gene_id "').rstrip('"'))
chess['Gene_name'] = chess['Gene_name'].map(lambda x: x.lstrip('gene_name "').rstrip('"'))

chess['Fasta_id'] = chess['Transcript_id'].map(str)+"|"+chess['Gene_id'].map(str)+"|"+chess['Gene_name'].map(str)+"|"+chess['Start'].map(str)+"|"+chess['End'].map(str)+'|'+chess['Type']


chess_transc = chess[chess['Type'] == 'transcript']

##make dict of transcript_id and fast id
keys = chess_transc['Transcript_id'].tolist()
values = chess_transc['Fasta_id'].tolist()
headers = dict(zip(keys, values))

from Bio import SeqIO

original_file = '/Users/mooreank/Desktop/Raph/requant/chess/chess.refseq.transcriptome.fa'
corrected_file = '/Users/mooreank/Desktop/Raph/requant/chess/chess.refseq.transcriptome.corrected.fa'

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        for key, value in headers.items():
            if key == record.id:
                record.id = value
                #print(record)
                SeqIO.write(record, corrected, 'fasta')
