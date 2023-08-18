import re
import os
import sys
import shutil
import zipfile
import requests
import threading
import pandas as pd

output_file_name = 'ncbi_input_species.csv'

# -------------------------READ API KEY --------------------------
api_key = str()
with open('ncbi_key.txt', 'r') as f:
    api_key = f.readline().strip()

if not api_key:
    print('Couldn\'t read API key in ncbi_key.txt')
    sys.exit(1)
# -----------------------------------------------------------------

if not os.path.exists('cds'):
    os.mkdir('cds')
if not os.path.exists('proteomes'):
    os.mkdir('proteomes')
if not os.path.exists('genomes'):
    os.mkdir('genomes')

ncbi_query_file = ''
with open('genome_result.txt', 'r') as f:
    ncbi_query_file = f.read()

organism_regex = r'^\d+\.\s(.*)\n'
matches = re.findall(organism_regex, ncbi_query_file, flags=re.MULTILINE)

def process_name(name):
    name = name.lower().replace(' ', '_')
    name = name.split('_')[:2]
    name = '_'.join(name)
    name = name.replace('[', '')
    name = name.replace(']', '')
    return name


def get_taxon_id(species_name, api_key) -> str:
    headers = {
        'Accept': 'application/json',
        'api-key': api_key,
    }

    response = requests.get(
        f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{species_name}',
        headers=headers
    )

    try:
        return response.json()['taxonomy_nodes'][0]['taxonomy']['tax_id']
    except KeyError:
        return ''


def get_accession(taxon_id, api_key) -> str:
    headers = {
        'Accept': 'application/json',
        'api-key': api_key,
    }

    response = requests.get(
        f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{taxon_id}/dataset_report?filters.reference_only=true&filters.assembly_source=all&filters.has_annotation=true&filters.exclude_paired_reports=true&filters.exclude_atypical=true&filters.assembly_version=current&filters.assembly_level=scaffold&filters.assembly_level=chromosome&filters.assembly_level=complete_genome&table_fields=assminfo-accession&table_fields=assminfo-name',
        headers=headers
    )

    try:
        return response.json()['reports'][0]['current_accession']
    except KeyError:
        return ''
    
# DOWNLOAD GENOME, PROTEOME, CDS
headers = {
    'Accept': 'application/zip',
}

params = {
    'include_annotation_type': [
        'GENOME_FASTA',
        'PROT_FASTA',
        'CDS_FASTA',
    ],
}

def fetch_url(name, names, genome_fasta_files, cds_fasta_files, original_names, allocated_i):
    new_name = process_name(name)
    file_name = name.lower().replace(' ', '_')
    file_name = file_name.replace('[', '')
    file_name = file_name.replace(']', '')

    print(new_name)

    cds_file_name = f'cds/{file_name}_cds.fna'
    protein_file_name = f'proteomes/{file_name}.faa'
    genome_file_name = f'genomes/{file_name}_genomic.fna'

    # TEST FOR EXISTENCE
    if not os.path.exists(cds_file_name) or not os.path.exists(genome_file_name) or not os.path.exists(protein_file_name):
        tax_id = get_taxon_id(species_name=name, api_key=api_key)  # GET TAXON
        accession = get_accession(taxon_id=tax_id, api_key=api_key)  # GET ACCESSION

        if tax_id == '':
            print(f'No taxonomy id found for {name}. Skipping...')
            return

        if accession == '':
            print(f'No accession id found for {name} taxon {tax_id}. Skipping...')
            return
        # ---------------------

        try:
            # DOWNLOAD
            response = requests.get(f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download', params=params, headers=headers)

            # WRITE DOWNLOAD
            with open(f'{name}{allocated_i}_download.zip', 'wb') as f:
                f.write(response.content)

            # EXTRACT DOWNLOADED FILE
            with zipfile.ZipFile(f'{name}{allocated_i}_download.zip', 'r') as zip_ref:
                os.mkdir(f'{name}{allocated_i}_download')
                zip_ref.extractall(f'{name}{allocated_i}_download')

            # NEW NAMES
            output_path = f'{name}{allocated_i}_download/ncbi_dataset/data/{accession}' # f'ncbi_dataset/data/{accession}'
            files = os.listdir(output_path)
            # --------------------------------------------------------------------

            # --------------MOVE FASTA FILES TO THE APPROPRIATE FOLDERS-----------
            for f in files:
                f_path = os.path.join(output_path, f)   # ncbi_dataset/data/{ACCESSION}/{f}

                if 'cds' in f:  # cds_from_genomic.fna
                    shutil.move(f_path, cds_file_name)  # cds/neurospora_crassa_cds.fna

                if 'protein' in f:  # protein.faa
                    shutil.move(f_path, protein_file_name)  # proteomes/neurospora_crassa.faa

                if accession in f:  # GCF_000182925.2_NC12 in GCF_000182925.2_NC12_genomic.fna
                    shutil.move(f_path, genome_file_name)  # genomes/neurospora_crassa_genomic.fna
            # --------------------------------------------------------------------
            # Clean up
            os.remove(f'{name}{allocated_i}_download.zip')
            shutil.rmtree(f'{name}{allocated_i}_download')

        except Exception as e:
            print(f'Something went wrong while downloading {name}')
            print(e)

            if os.path.exists(f'{name}{allocated_i}_download.zip'):
                os.remove(f'{name}{allocated_i}_download.zip')

            if os.path.exists(cds_file_name):
                os.remove(cds_file_name)
            
            if os.path.exists(protein_file_name):
                os.remove(protein_file_name)

            if os.path.exists(genome_file_name):
                os.remove(genome_file_name)

            if os.path.exists(f'{name}{allocated_i}'):
                shutil.rmtree(f'{name}{allocated_i}')

            return

    names[allocated_i] = new_name
    genome_fasta_files[allocated_i] = f'{file_name}_genomic.fna'
    cds_fasta_files[allocated_i] = f'{file_name}_cds.fna'
    original_names[allocated_i] = file_name


def fetch_url_chunk(chunk):
    threads = []

    names = [None] * len(chunk)
    genome_fasta_files = [None] * len(chunk)
    cds_fasta_files = [None] * len(chunk)
    original_names = [None] * len(chunk)

    for i, name in enumerate(chunk):
        thread = threading.Thread(
            target=fetch_url,
            args=(name, names, genome_fasta_files, cds_fasta_files, original_names, i)
            )
        
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

chunk_size = 15
chunks = [matches[i:i+chunk_size] for i in range(0, len(matches), chunk_size)]

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