import os
import re
import gzip
import shutil
import requests
import threading
import pandas as pd

df = pd.read_csv('ensemblFungi.csv')
df = df.drop(columns=df.columns[0])

if not os.path.exists('genomes'):
    os.mkdir('genomes')
if not os.path.exists('cds'):
    os.mkdir('cds')
if not os.path.exists('proteomes'):
    os.mkdir('proteomes')

output_file_name = 'ensembl_input_species.csv'

def process_name(name):
    name = name.lower().replace(' ', '_')
    name = name.split('_')[:2]
    name = '_'.join(name)
    name = name.replace('[', '')
    name = name.replace(']', '')
    return name


def fetch_url(row, names, genome_fasta_files, cds_fasta_files, original_names, allocated_i):
    name = row['Species']
    name = name.lower().replace(' ', '_')
    name = name.replace('[', '')
    name = name.replace(']', '')

    new_name = process_name(name)
    
    print(new_name)

    cds_file_name = f'cds/{name}_cds.fna'
    protein_file_name = f'proteomes/{name}.faa'
    genome_file_name = f'genomes/{name}_genomic.fna'

    cds_url = row['cds_url']
    dna_url = row['dna_url']
    prot_url = row['prot_url']

    # -- CDS --
    if not os.path.exists(cds_file_name) or not os.path.exists(genome_file_name) or not os.path.exists(protein_file_name):
        try:
            cds_r = requests.get(cds_url)
            if str(cds_r.content[0:15]) == "b'<!doctype html>'":
                print(f'{name} no cds')
                return

            m = re.findall(r"<a href=\"([^?/]+)\"", str(cds_r.content))
            index = [idx for idx, s in enumerate(m) if '.gz' in s][0]
            file = requests.get(cds_url + m[index])

            open(f'{name}_cds.gz', 'wb').write(file.content)
            with gzip.open(f'{name}_cds.gz', 'rb') as f_in:
                with open(cds_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f'{name}_cds.gz')
            
            # -- DNA --
            dna_r = requests.get(dna_url)
            if str(dna_r.content[0:15]) == "b'<!doctype html>'":
                print('no genome')
                os.remove(f'{name}_cds.gz')
                return

            m = re.findall(r"<a href=\"([^?/]+dna.toplevel.fa.gz)\"", str(dna_r.content))
            index = [idx for idx, s in enumerate(m) if '.gz' in s][0]
            file = requests.get(dna_url + m[index])

            open(f'{name}_dna.gz', 'wb').write(file.content)
            with gzip.open(f'{name}_dna.gz', 'rb') as f_in:
                with open(genome_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f'{name}_dna.gz')
        
            # -- PROT --
            prot_r = requests.get(prot_url)
            if str(prot_r.content[0:15]) == "b'<!doctype html>'":
                print(f'{name} no prot')
                os.remove(f'{name}_cds.gz')
                os.remove(f'{name}_dna.gz')
                return
            
            m = re.findall(r"<a href=\"([^?/]+pep.all.fa.gz)\"", str(prot_r.content))
            index = [idx for idx, s in enumerate(m) if '.gz' in s][0]
            file = requests.get(prot_url + m[index])

            open(f'{name}_prot.gz', 'wb').write(file.content)
            with gzip.open(f'{name}_prot.gz', 'rb') as f_in:
                with open(protein_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f'{name}_prot.gz')
        
        except Exception as e:
            print(f'Something went wrong while downloading {name}')
            print(e)

            if os.path.exists(f'{cds_file_name}.gz'):
                os.remove(f'{cds_file_name}.gz')

            if os.path.exists(f'{genome_file_name}.gz'):
                os.remove(f'{genome_file_name}.gz')

            if os.path.exists(f'{protein_file_name}.gz'):
                os.remove(f'{protein_file_name}.gz')

            if os.path.exists(f'{name}_dna.gz'):
                os.remove(f'{name}_dna.gz')

            if os.path.exists(f'{name}_cds.gz'):
                os.remove(f'{name}_cds.gz')

            if os.path.exists(f'{name}_prot.gz'):
                os.remove(f'{name}_prot.gz')

            names[allocated_i] = None
            genome_fasta_files[allocated_i] = None
            cds_fasta_files[allocated_i] = None
            original_names[allocated_i] = None
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

chunk_size = 10
chunks = [df.iloc[i:i+chunk_size] for i in range(0, len(df), chunk_size)]

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