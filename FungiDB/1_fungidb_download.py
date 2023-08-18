import os
import requests
import threading
import pandas as pd

fungidb = pd.read_csv('fungidb.csv', sep='\t')
fungidb = fungidb[fungidb['Is Reference Strain'] == 'yes'].reset_index(drop=True)

fungidb['short_name'] = fungidb.Species.str.replace(' ', '_').str.lower().str.split('_').str[:2].str.join('_').str.replace('[', '', regex=True).replace(']', '', regex=True)

# fungidb = fungidb.drop_duplicates(subset='short_name').reset_index(drop=True)

manifest = pd.read_csv('GenomeDataTypes_Summary.csv')
manifest['short_name'] = manifest['Species'].str.replace(' ', '_').str.lower().str.split('_').str[:2].str.join('_').str.replace('[', '', regex=True).replace(']', '', regex=True)

shared = manifest[manifest['Species'].isin(fungidb['Species'])]
shared = shared[shared['Is Reference Strain'] == 'yes']
shared = shared[shared['Protein coding genes'] != 0]
# shared = shared[~shared.duplicated(subset='short_name')].reset_index(drop=True)

if not os.path.exists('genomes'):
    os.mkdir('genomes')
if not os.path.exists('cds'):
    os.mkdir('cds')
if not os.path.exists('proteomes'):
    os.mkdir('proteomes')

output_file_name = 'fungidb_input_species.csv'

def fetch_url(row, names, genome_fasta_files, cds_fasta_files, original_names, allocated_i):
    new_name = row['short_name']
    name = row['Species']
    name = name.lower().replace(' ', '_')
    name = name.replace('[', '')
    name = name.replace(']', '')

    print(new_name)
    
    cds_file_name = f'cds/{name}_cds.fna'
    protein_file_name = f'proteomes/{name}.faa'
    genome_file_name = f'genomes/{name}_genomic.fna'

    protein_link = row['Protein Fasta Download Link']
    genome_link = row['Genome Fasta Download Link']
    cds_link = protein_link.replace('_AnnotatedProteins.fasta', '_AnnotatedCDSs.fasta')

    if not os.path.exists(cds_file_name) or not os.path.exists(protein_file_name) or not os.path.exists(genome_file_name):
        try:
            cds_r = requests.get(cds_link)
            if str(cds_r.content[0:15]) == "b'<!doctype html>'":
                print('no cds')
                return
            
            open(cds_file_name, 'wb').write(cds_r.content)

            prot_r = requests.get(protein_link)
            if str(prot_r.content[0:15]) == "b'<!doctype html>'":
                print('no prot')
                os.remove(cds_file_name)
                return
            
            open(protein_file_name, 'wb').write(prot_r.content)

            genome_r = requests.get(genome_link)
            if str(genome_r.content[0:15]) == "b'<!doctype html>'":
                print('no genome')
                os.remove(cds_file_name)
                os.remove(protein_file_name)
                return
            
            open(genome_file_name, 'wb').write(genome_r.content)

        except Exception as e:
            print(f'Something went wrong while downloading {name}')
            print(e)

            if os.path.exists(f'{cds_file_name}.gz'):
                os.remove(f'{cds_file_name}.gz')

            if os.path.exists(f'{genome_file_name}.gz'):
                os.remove(f'{genome_file_name}.gz')

            if os.path.exists(f'{protein_file_name}.gz'):
                os.remove(f'{protein_file_name}.gz')

            return

    names[allocated_i] = new_name
    genome_fasta_files[allocated_i] = f'{name}_genomic.fna'
    cds_fasta_files[allocated_i] = f'{name}_cds.fna'
    original_names[allocated_i] = name


def fetch_url_chunk(rows):
    threads = []

    names = [None] * len(rows)
    genome_fasta_files = [None] * len(rows)
    cds_fasta_files = [None] * len(rows)
    original_names = [None] * len(rows)

    i = 0
    for _, row in rows.iterrows():
        thread = threading.Thread(
            target=fetch_url,
            args=(row, names, genome_fasta_files, cds_fasta_files, original_names, i)
            )
        
        i += 1
        
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    names = [name for name in names if name is not None]
    genome_fasta_files = [name for name in genome_fasta_files if name is not None]
    cds_fasta_files = [name for name in cds_fasta_files if name is not None]
    original_names = [name for name in original_names if name is not None]

    return names, genome_fasta_files, cds_fasta_files, original_names


all_names = list()
all_genome_fasta_files = list()
all_cds_fasta_files = list()
all_original_names = list()

chunk_size = 5
chunks = [shared.iloc[i:i+chunk_size] for i in range(0, len(shared), chunk_size)]

for url_chunk in chunks:
    names, genome_fasta_files, cds_fasta_files, original_names = fetch_url_chunk(url_chunk)
    all_names.extend(names)
    all_genome_fasta_files.extend(genome_fasta_files)
    all_cds_fasta_files.extend(cds_fasta_files)
    all_original_names.extend(original_names)


with open('names.txt', 'w') as f:
    for n in all_names:
        f.writelines(n + '\n')

df = pd.DataFrame({
    'species_name': all_names,
    'genome_file_name': all_genome_fasta_files,
    'cds_file_name': all_cds_fasta_files,
    'original_name': all_original_names
})

df.to_csv(output_file_name, index=False)