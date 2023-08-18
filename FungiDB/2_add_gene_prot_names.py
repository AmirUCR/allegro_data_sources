import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

gene_prot_regex = r'transcript=([^\s]*)'

i = 0
for cds_f in os.listdir('cds'):
    cds_records = list()
    cds_path = os.path.join('cds', cds_f)
    records = list(SeqIO.parse(open(cds_path), 'fasta'))

    for record in records:
        new_record = SeqRecord(
            seq=record.seq,
            id=record.id,
            description=(record.description + f' | [protein_id={record.id}]')
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
        prot_match = re.search(gene_prot_regex, record.description)

        if prot_match:
            new_record = SeqRecord(
                seq=record.seq,
                id=prot_match.group(1),
                description=str(record.id) + ' | ' + str(record.description)
                )
            prot_records.append(new_record)
    
    with open(prot_path, 'w') as f:
        SeqIO.write(prot_records, f, 'fasta')

    if i % 20 == 0: print('[Protein] Done with', i, 'species.', end='\r')
    i += 1