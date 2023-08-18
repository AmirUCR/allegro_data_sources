import os
import shutil
import pandas as pd

fdb_dir = 'FungiDB'
ncbi_dir = 'NCBI'
# onek_dir = '1000_Fungi_Project'
ensembl_dir = 'EnsemblFungi'
mycocosm_dir = 'Mycocosm'


fdb = pd.read_csv(f'{fdb_dir}/fungidb_input_species.csv')
fdb['source'] = 'fdb'

ncbi = pd.read_csv(f'{ncbi_dir}/ncbi_input_species.csv')
ncbi['source'] = 'ncbi'

# onek = pd.read_csv(f'{onek_dir}/jgi_input_species.csv')
# onek['source'] = 'onek'
mycocosm = pd.read_csv(f'{mycocosm_dir}/mycocosm_input_species.csv')
mycocosm['source'] = 'mycocosm'

ensembl = pd.read_csv(f'{ensembl_dir}/ensembl_input_species.csv')
ensembl['source'] = 'ensembl'

fdb = fdb.drop_duplicates(subset='species_name').reset_index(drop=True)
ncbi = ncbi.drop_duplicates(subset='species_name').reset_index(drop=True)
# onek = onek.drop_duplicates(subset='species_name').reset_index(drop=True)
ensembl = ensembl.drop_duplicates(subset='species_name').reset_index(drop=True)
mycocosm = mycocosm.drop_duplicates(subset='species_name').reset_index(drop=True)

concat = pd.concat([ncbi, fdb, ensembl, mycocosm], ignore_index=True)
all_dupes = concat[concat.duplicated(subset=['species_name'])]

concat = concat[~concat.duplicated(subset=['species_name'])].reset_index(drop=True)

concat_destination_dir = 'fourdbs_concat'

if not os.path.exists(concat_destination_dir):
    os.mkdir(concat_destination_dir)

if not os.path.exists(os.path.join(concat_destination_dir, 'cds')):
    os.makedirs(os.path.join(concat_destination_dir, 'cds'))

if not os.path.exists(os.path.join(concat_destination_dir, 'genomes')):
    os.makedirs(os.path.join(concat_destination_dir, 'genomes'))

if not os.path.exists(os.path.join(concat_destination_dir, 'proteomes')):
    os.makedirs(os.path.join(concat_destination_dir, 'proteomes'))

for idx, row in concat.iterrows():
    cds_f_name = row['cds_file_name']
    genome_f_name = row['genome_file_name']
    protein_f_name = row['original_name'] + '.faa'
    
    if row['source'] == 'mycocosm':
        shutil.copy(
            os.path.join(mycocosm_dir, 'cds', cds_f_name),
            os.path.join(concat_destination_dir, 'cds'))
        shutil.copy(
            os.path.join(mycocosm_dir, 'genomes', genome_f_name),
            os.path.join(concat_destination_dir, 'genomes'))
        shutil.copy(
            os.path.join(mycocosm_dir, 'proteomes', protein_f_name),
            os.path.join(concat_destination_dir, 'proteomes'))

    if row['source'] == 'fdb':
        shutil.copy(
            os.path.join(fdb_dir, 'cds', cds_f_name),
            os.path.join(concat_destination_dir, 'cds'))
        shutil.copy(
            os.path.join(fdb_dir, 'genomes', genome_f_name),
            os.path.join(concat_destination_dir, 'genomes'))
        shutil.copy(
            os.path.join(fdb_dir, 'proteomes', protein_f_name),
            os.path.join(concat_destination_dir, 'proteomes'))

    if row['source'] == 'ensembl':
        shutil.copy(
            os.path.join(ensembl_dir, 'cds', cds_f_name),
            os.path.join(concat_destination_dir, 'cds'))
        shutil.copy(
            os.path.join(ensembl_dir, 'genomes', genome_f_name),
            os.path.join(concat_destination_dir, 'genomes'))
        shutil.copy(
            os.path.join(ensembl_dir, 'proteomes', protein_f_name),
            os.path.join(concat_destination_dir, 'proteomes'))

    elif row['source'] == 'ncbi':
        shutil.copy(
            os.path.join(ncbi_dir, 'cds', cds_f_name),
            os.path.join(concat_destination_dir, 'cds'))
        shutil.copy(
            os.path.join(ncbi_dir, 'genomes', genome_f_name),
            os.path.join(concat_destination_dir, 'genomes'))
        shutil.copy(
            os.path.join(ncbi_dir, 'proteomes', protein_f_name),
            os.path.join(concat_destination_dir, 'proteomes'))

concat.to_csv(
    os.path.join(concat_destination_dir, 'fourdbs_input_species.csv'),
    index=False)