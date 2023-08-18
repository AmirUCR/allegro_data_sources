import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

gene_prot_regex = r'[^|]+$'

i = 0
for cds_f in os.listdir('cds'):
    cds_records = list()
    cds_path = os.path.join('cds', cds_f)
    records = list(SeqIO.parse(open(cds_path), 'fasta'))

    for record in records:
        gene_match = re.search(gene_prot_regex, record.id)

        if gene_match:
            new_record = SeqRecord(
                seq=record.seq,
                id=record.id,
                description=(record.description + ' [protein_id=' + gene_match.group(0) + ']')
                )
            cds_records.append(new_record)
    
    with open(cds_path, 'w') as f:
        SeqIO.write(cds_records, f, 'fasta')
    
    if i % 20 == 0: print('[CDS] Done with', i, 'species.', end='\r')
    i += 1

i = 0
for prot_f in os.listdir('proteomes'):
    prot_records = list()
    prot_path = os.path.join('proteomes', prot_f)
    records = list(SeqIO.parse(open(prot_path), 'fasta'))

    for record in records:
        prot_match = re.search(gene_prot_regex, record.id)

        if prot_match:
            new_record = SeqRecord(
                seq=record.seq,
                id=prot_match.group(0),
                description=str(record.id)
                )
            prot_records.append(new_record)
    
    with open(prot_path, 'w') as f:
        SeqIO.write(prot_records, f, 'fasta')

    if i % 20 == 0: print('[Protein] Done with', i, 'species.', end='\r')
    i += 1