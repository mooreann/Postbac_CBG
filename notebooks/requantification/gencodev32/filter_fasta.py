#!/usr/bin/env python

#import sys
#from sets import Set
from Bio import SeqIO

WRKDIR = '/home/mooreank/salmon/'
ids = '{}bad.gencode.v32.transcript_ids.txt'.format(WRKDIR)
#original_fasta = '{}gencode.v32.transcripts.fa'.format(WRKDIR)
#corrected_fasta = '{}filtered.gencode.v32.transcripts.fa'.format(WRKDIR)

# read the first file given and generate a set (faster iteration respect lists
##get list of transcript ids to remove
identifiers = []

with open(ids, 'r') as fi:
    for line in fi:
        line = line.strip()
        identifiers.append(str(line).replace(">", ""))
        
##get list of all record ids in original fasta
record_ids = []

with open('{}gencode.v32.transcripts.fa'.format(WRKDIR)) as original_fasta:
    records = SeqIO.parse(original_fasta, 'fasta')
    for record in records:
        record_ids.append(record.id)



##get list of record ids to keep
filtered = [i for i in record_ids if not any(i for j in identifiers if str(j) in i)]


##run through original fasta and keep only records with record ids in filtered list
with open('{}gencode.v32.transcripts.fa'.format(WRKDIR)) as original_fasta, open('{}/nohup.filtered.gencode.v32.transcripts.fa'.format(WRKDIR), 'w') as corrected_fasta:

        records = SeqIO.parse(original_fasta, 'fasta')
        for record in records:
            for x in filtered:
                if x == record.id:
                    SeqIO.write(record, corrected_fasta, 'fasta')

