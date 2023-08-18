import os
import gzip
import shutil
import requests
import threading
import pandas as pd
import xml.etree.ElementTree as ET
from http.cookiejar import MozillaCookieJar

# ----- REPLACE THESE WITH YOUR JGI USERNAME AND PASSWORD --------
data = {
    'login': 'YOUR_JGI_USERNAME',  # Your JGI account username
    'password': 'YOUR_JGI_PASSWORD',  # Your JGI account password
}
# ----- You do not need to modify anything after this point ------

if data['login'] == 'YOUR_JGI_USERNAME' and data['password'] == 'YOUR_JGI_PASSWORD':
    with open('authenticate.txt', 'r') as f:
        data['login'] = f.readline().strip()
        data['password'] = f.readline().strip()

# Clean up previous run.
if os.path.exists('cookies'):
    os.remove('cookies')

cookies = MozillaCookieJar()
with requests.Session() as session:
    session.cookies = cookies
    response = session.post('https://signon.jgi.doe.gov/signon/create', data=data)
    cookies.save('cookies', ignore_discard=True, ignore_expires=True)

# Create directories if nonexistent.
if not os.path.exists('genomes'):
    os.mkdir('genomes')
if not os.path.exists('cds'):
    os.mkdir('cds')
if not os.path.exists('proteomes'):
    os.mkdir('proteomes')

# Name of the output file.
output_file_name = 'mycocosm_input_species.csv'

def process_name(name):
    name = name.lower().replace(' ', '_')
    name = name.split('_')[:2]
    name = '_'.join(name)
    name = name.replace('[', '')
    name = name.replace(']', '')
    return name

# Read the file tree
tree = ET.parse('get-directory.xml')
root = tree.getroot()

assembly_children = root.find(
'''.//folder[@name='Assembly']/folder[@name='Genome Assembly (masked)']'''
)

annotated_cds_children = [child for child in root.find(
'''.//folder[@name='Annotation']/folder[@name='Filtered Models (\"best\")']/folder[@name='CDS']'''
) if 'GeneCatalog' in child.attrib['url'] ]

annotated_protein_children = [child for child in root.find(
'''.//folder[@name='Annotation']/folder[@name='Filtered Models (\"best\")']/folder[@name='Proteins']'''
) if 'GeneCatalog' in child.attrib['url'] ]

class Fungus:
    def __init__(self, name='', cds_url='', protein_url='', genome_url=''):
        self.name = name
        self.cds_url = cds_url
        self.protein_url = protein_url
        self.genome_url = genome_url


fungi_dict = dict()

if assembly_children:
    for idx, child in enumerate(assembly_children):
        f = None
        name = child.get('label')
        if name:

            # Only keep the first occurance
            if name in fungi_dict:
                continue

            url = child.get('url')
            if url:
                f = Fungus()
                f.name = name
                f.genome_url = url
                fungi_dict[name] = f
            else:
                print('No assembly URL attribute found for', idx, '. Skipping.')
                continue
        else:
            print('No label found for assembly entry', idx, '.')


if annotated_cds_children:
    for idx, child in enumerate(annotated_cds_children):
        name = child.get('label')

        if name and name in fungi_dict:
                url = child.get('url')
                
                if url:
                    fungi_dict[name].cds_url = url
                else:
                    print('No CDS URL attribute found for', idx, '. Skipping.')
                    del fungi_dict[name]
                    continue
        else:
            print('No label attribute found for CDS entry', idx, '.')
            continue


if annotated_protein_children:
    for idx, child in enumerate(annotated_protein_children):
        name = child.get('label')

        if name and name in fungi_dict:
                url = child.get('url')
                
                if url:
                    fungi_dict[name].protein_url = url
                else:
                    print('No prot URL attribute found for', idx, '. Skipping.')
                    del fungi_dict[name]
                    continue
        else:
            print('No label attribute found for protein entry', idx, '.')
            continue

def fetch_url(name, cds_url, genome_url, protein_url,names, genome_fasta_files, cds_fasta_files, original_names, allocated_i):
    new_name = process_name(name)
    name = name.lower().replace(' ', '_')
    name = name.replace('[', '')
    name = name.replace(']', '')

    print(new_name)

    cds_file_name = f'cds/{name}_cds.fna'
    protein_file_name = f'proteomes/{name}.faa'
    genome_file_name = f'genomes/{name}_genomic.fna'

    if not os.path.exists(cds_file_name) or not os.path.exists(protein_file_name) or not os.path.exists(genome_file_name):
        try:
            # GET CDS
            response = requests.get(
                'https://genome.jgi.doe.gov{URL}'.format(URL=cds_url),
                cookies=cookies,
            )

            with open(f'{cds_file_name}.gz', 'wb') as f:
                f.write(response.content)

            with gzip.open(f'{cds_file_name}.gz', 'rb') as f_in:
                with open(cds_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            os.remove(f'{cds_file_name}.gz')

            # GET PROT
            response = requests.get(
                'https://genome.jgi.doe.gov{URL}'.format(URL=protein_url),
                cookies=cookies,
            )

            with open(f'{protein_file_name}.gz', 'wb') as f:
                f.write(response.content)

            with gzip.open(f'{protein_file_name}.gz', 'rb') as f_in:
                with open(protein_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            os.remove(f'{protein_file_name}.gz')

            # GET DNA
            response = requests.get(
                'https://genome.jgi.doe.gov{URL}'.format(URL=genome_url),
                cookies=cookies,
            )

            with open(f'{genome_file_name}.gz', 'wb') as f:
                f.write(response.content)

            with gzip.open(f'{genome_file_name}.gz', 'rb') as f_in:
                with open(genome_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            os.remove(f'{genome_file_name}.gz')
            
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
    

def fetch_url_chunk(chunk):
    threads = []

    names = [None] * len(chunk)
    genome_fasta_files = [None] * len(chunk)
    cds_fasta_files = [None] * len(chunk)
    original_names = [None] * len(chunk)

    for i, fungi_obj in enumerate(chunk):
        name = fungi_obj.name
        cds_url = fungi_obj.cds_url
        genome_url = fungi_obj.genome_url
        protein_url = fungi_obj.protein_url

        thread = threading.Thread(
            target=fetch_url,
            args=(name, cds_url, genome_url, protein_url, names, genome_fasta_files, cds_fasta_files, original_names, i)
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

chunk_size = 1
chunks = [list(fungi_dict.values())[i:i+chunk_size] for i in range(0, len(fungi_dict), chunk_size)]

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