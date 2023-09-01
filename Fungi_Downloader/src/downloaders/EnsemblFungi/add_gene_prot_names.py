import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def fix_ids():
    i = 0
    for cds_f in os.listdir('data/EnsemblFungi/cds'):
        cds_records = list()
        cds_path = os.path.join('data/EnsemblFungi/cds', cds_f)
        records = list(SeqIO.parse(open(cds_path), 'fasta'))

        for record in records:
            new_record = SeqRecord(
                seq=record.seq,
                id=record.id,
                description=(record.description + f' [protein_id={record.id}]')
                )
            cds_records.append(new_record)
        
        with open(cds_path, 'w') as f:
            SeqIO.write(cds_records, f, 'fasta')
        
        if i % 20 == 0: print('[CDS] Done with', i, 'species.', end='\r')
        i += 1
    print()