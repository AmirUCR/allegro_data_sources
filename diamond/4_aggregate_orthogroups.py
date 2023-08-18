import os
import pandas

outputs_directory = 'outputs'

tsv_path = os.path.join(outputs_directory, 'saccharomyces_cerevisiae.tsv')

headers = ['query_accession', 'target_accession', 'sequence_identity', 'length', 'mismatches',
           'gap_openings', 'query_start', 'query_end', 'target_start', 'target_end', 'e_value',
           'bit_score']

ref_df = pandas.read_csv(tsv_path, names=headers, sep='\t')

ref_gene_accessions = ref_df['query_accession']

non_ref_species = sorted([x.split('.tsv')[0] for x in os.listdir(outputs_directory)])

orthogroups_dict = {
    'Orthogroup': ref_gene_accessions.values
    }

for nrs in non_ref_species:
    orthogroups_dict[nrs] = list()

    nrs_path = os.path.join(outputs_directory, nrs + '.tsv')

    nrs_df = pandas.read_csv(nrs_path, names=headers, sep='\t')

    for ref_gene_accession in ref_gene_accessions:
        target_accession = 'Not found'
        if ref_gene_accession in nrs_df['query_accession'].values:
            target_accession = nrs_df[nrs_df['query_accession'] == ref_gene_accession]['target_accession'].values[0]

        orthogroups_dict[nrs].append(target_accession)

pandas.DataFrame.from_dict(orthogroups_dict).to_csv('Orthogroups.tsv', sep='\t', index=False)

print('Done. Check Orthogroups.tsv')