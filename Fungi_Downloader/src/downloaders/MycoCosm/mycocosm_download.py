import os
import gzip
import shutil
import requests
import threading
import pandas as pd
import xml.etree.ElementTree as ET
from http.cookiejar import MozillaCookieJar

from downloaders.MycoCosm.add_gene_prot_names import fix_ids
from utils.path_generator import generate_dirs
from utils.name_processor import process_two_part_name, process_name


# Helper class
class Fungus:
    def __init__(self, name='', cds_url='', protein_url='', genome_url=''):
        self.name = name
        self.cds_url = cds_url
        self.protein_url = protein_url
        self.genome_url = genome_url


class MycoCosm_Downloader:
    def __init__(self) -> None:

        # ----- REPLACE THESE WITH YOUR JGI USERNAME AND PASSWORD --------
        data = {
            'login': 'YOUR_JGI_USERNAME',  # Your JGI account username
            'password': 'YOUR_JGI_PASSWORD',  # Your JGI account password
        }
        # ----- You do not need to modify anything after this point ------

        if data['login'] == 'YOUR_JGI_USERNAME' and data['password'] == 'YOUR_JGI_PASSWORD':
            with open('src/downloaders/MycoCosm/authenticate.txt', 'r') as f:
                data['login'] = f.readline().strip()
                data['password'] = f.readline().strip()

        # Clean up previous run.
        if os.path.exists('src/downloaders/MycoCosm/cookies'):
            os.remove('src/downloaders/MycoCosm/cookies')

        self.cookies = MozillaCookieJar()
        with requests.Session() as session:
            session.cookies = self.cookies
            response = session.post('https://signon.jgi.doe.gov/signon/create', data=data)
            self.cookies.save('src/downloaders/MycoCosm/cookies', ignore_discard=True, ignore_expires=True)

        generate_dirs('data/MycoCosm')

        # Name of the output file.
        self.output_file_name = 'data/MycoCosm/mycocosm_input_species.csv'

        # Read the file tree
        tree = ET.parse('src/downloaders/MycoCosm/get-directory.xml')
        root = tree.getroot()

        assembly_children = root.find(
        '''.//folder[@name='Assembly']/folder[@name='Genome Assembly (masked)']'''
        )

        annotated_cds_children = [child for child in root.find(
        '''.//folder[@name='Annotation']/folder[@name='Filtered Models (\"best\")']/folder[@name='CDS']'''
        ) if 'GeneCatalog' in child.attrib['url'] ]

        annotated_protein_children = [child for child in root.find(
        '''.//folder[@name='Annotation']/folder[@name='Filtered Models (\"best\")']/folder[@name='Proteins']'''
        ) if 'GeneCatalog' in child.attrib['url'] and '.aa.fasta.gz' in child.attrib['url'] ]

        self.fungi_dict = dict()

        if assembly_children:
            for idx, child in enumerate(assembly_children):
                f = None
                name = child.get('label')
                if name:

                    # Only keep the first occurance
                    if name in self.fungi_dict:
                        continue

                    url = child.get('url')
                    if len(url) > 10:
                        f = Fungus()
                        f.name = name
                        f.genome_url = url
                        self.fungi_dict[name] = f
                    else:
                        print('No assembly URL attribute found for', idx, '. Skipping.')
                        continue
                else:
                    print('No label found for assembly entry', idx, '.')


        if annotated_cds_children:
            for idx, child in enumerate(annotated_cds_children):
                name = child.get('label')

                if name and name in self.fungi_dict:
                    url = child.get('url')
                    
                    if len(url) > 10:
                        self.fungi_dict[name].cds_url = url
                    else:
                        print('No CDS URL attribute found for', idx, '. Skipping.')
                        del self.fungi_dict[name]
                        continue
                else:
                    print('No label attribute found for CDS entry', idx, '.')
                    continue


        if annotated_protein_children:
            for idx, child in enumerate(annotated_protein_children):
                name = child.get('label')

                if name and name in self.fungi_dict:
                    url = child.get('url')
                    
                    if len(url) > 10:
                        self.fungi_dict[name].protein_url = url
                    else:
                        print('No prot URL attribute found for', idx, '. Skipping.')
                        del self.fungi_dict[name]
                        continue
                else:
                    print('No label attribute found for protein entry', idx, '.')
                    continue

        # Remove entries with missing URLs
        for name, obj in list(self.fungi_dict.items()):
            if len(obj.cds_url) < 1 or len(obj.genome_url) < 1 or len(obj.protein_url) < 1:
                del self.fungi_dict[name]


    def fetch_url(self, name, cds_url, genome_url, protein_url, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, allocated_i):
        new_name = process_two_part_name(name)
        file_name = process_name(name)

        print(new_name)

        cds_file_name = f'data/MycoCosm/cds/{file_name}_cds.fna'
        protein_file_name = f'data/MycoCosm/proteomes/{file_name}.faa'
        genome_file_name = f'data/MycoCosm/genomes/{file_name}_genomic.fna'

        if not os.path.exists(cds_file_name) or not os.path.exists(protein_file_name) or not os.path.exists(genome_file_name):
            try:
                # GET CDS
                response = requests.get(f'https://genome.jgi.doe.gov{cds_url}', cookies=self.cookies)

                with open(f'{cds_file_name}_mycocosm.gz', 'wb') as f:
                    f.write(response.content)

                with gzip.open(f'{cds_file_name}_mycocosm.gz', 'rb') as f_in:
                    with open(cds_file_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(f'{cds_file_name}_mycocosm.gz')

                # GET PROT
                response = requests.get(f'https://genome.jgi.doe.gov{protein_url}', cookies=self.cookies)

                with open(f'{protein_file_name}_mycocosm.gz', 'wb') as f:
                    f.write(response.content)

                with gzip.open(f'{protein_file_name}_mycocosm.gz', 'rb') as f_in:
                    with open(protein_file_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(f'{protein_file_name}_mycocosm.gz')

                # GET DNA
                response = requests.get(f'https://genome.jgi.doe.gov{genome_url}', cookies=self.cookies)

                with open(f'{genome_file_name}_mycocosm.gz', 'wb') as f:
                    f.write(response.content)

                with gzip.open(f'{genome_file_name}_mycocosm.gz', 'rb') as f_in:
                    with open(genome_file_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(f'{genome_file_name}_mycocosm.gz')
                
            except Exception as e:
                print(f'Something went wrong while downloading {name}')
                print(e)

                if os.path.exists(f'{cds_file_name}_mycocosm.gz'):
                    os.remove(f'{cds_file_name}_mycocosm.gz')

                if os.path.exists(f'{genome_file_name}_mycocosm.gz'):
                    os.remove(f'{genome_file_name}_mycocosm.gz')

                if os.path.exists(f'{protein_file_name}_mycocosm.gz'):
                    os.remove(f'{protein_file_name}_mycocosm.gz')

                return

        names[allocated_i] = new_name
        genome_fasta_files[allocated_i] = f'{file_name}_genomic.fna'
        cds_fasta_files[allocated_i] = f'{file_name}_cds.fna'
        original_names[allocated_i] = file_name
        cds_urls[allocated_i] = f"curl -c cookies.txt -d \"login=YOUR_JGI_USERNAME&password=YOUR_JGI_PASSWORD\" \"https://signon.jgi.doe.gov/signon/create && curl -b cookies.txt \"https://genome.jgi.doe.gov{cds_url}\" --output {cds_file_name}_cds_download.gz"
        genome_urls[allocated_i] = f"curl -c cookies.txt -d \"login=YOUR_JGI_USERNAME&password=YOUR_JGI_PASSWORD\" \"https://signon.jgi.doe.gov/signon/create && curl -b cookies.txt \"https://genome.jgi.doe.gov{genome_url}\" --output {genome_file_name}_genome_download.gz"
        proteome_urls[allocated_i] = f"curl -c cookies.txt -d \"login=YOUR_JGI_USERNAME&password=YOUR_JGI_PASSWORD\" \"https://signon.jgi.doe.gov/signon/create && curl -b cookies.txt \"https://genome.jgi.doe.gov{protein_url}\" --output {protein_file_name}_prot_download.gz"


    def fetch_url_chunk(self, chunk):
        threads = []

        names = [None] * len(chunk)
        genome_fasta_files = [None] * len(chunk)
        cds_fasta_files = [None] * len(chunk)
        original_names = [None] * len(chunk)
        cds_urls = [None] * len(chunk)
        genome_urls = [None] * len(chunk)
        proteome_urls = [None] * len(chunk)

        for i, fungi_obj in enumerate(chunk):
            name = fungi_obj.name
            cds_url = fungi_obj.cds_url
            genome_url = fungi_obj.genome_url
            protein_url = fungi_obj.protein_url

            thread = threading.Thread(
                target=self.fetch_url,
                args=(name, cds_url, genome_url, protein_url, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, i)
                )
            
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

        names = [name for name in names if name is not None]
        genome_fasta_files = [name for name in genome_fasta_files if name is not None]
        cds_fasta_files = [name for name in cds_fasta_files if name is not None]
        original_names = [name for name in original_names if name is not None]
        cds_urls = [url for url in cds_urls if url is not None]
        genome_urls = [url for url in genome_urls if url is not None]
        proteome_urls = [url for url in proteome_urls if url is not None]

        return names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls


    def download(self):
        all_names = list()
        all_genome_fasta_files = list()
        all_cds_fasta_files = list()
        all_original_names = list()
        all_cds_urls = list()
        all_genome_urls = list()
        all_proteome_urls = list()

        chunk_size = 1
        chunks = [list(self.fungi_dict.values())[i:i+chunk_size] for i in range(0, len(self.fungi_dict), chunk_size)]

        for url_chunk in chunks:
            names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls = self.fetch_url_chunk(url_chunk)
            all_names.extend(names)
            all_genome_fasta_files.extend(genome_fasta_files)
            all_cds_fasta_files.extend(cds_fasta_files)
            all_original_names.extend(original_names)
            all_cds_urls.extend(cds_urls)
            all_genome_urls.extend(genome_urls)
            all_proteome_urls.extend(proteome_urls)

        with open('data/MycoCosm/names.txt', 'w') as f:
            for n in all_names:
                f.writelines(n + '\n')

        df = pd.DataFrame({
            'species_name': all_names,
            'genome_file_name': all_genome_fasta_files,
            'cds_file_name': all_cds_fasta_files,
            'original_name': all_original_names,
            'cds_url': all_cds_urls,    
            'genome_url': all_genome_urls,
            'proteome_url': all_proteome_urls
        })

        df.to_csv(self.output_file_name, index=False)
        
        print('Fixing Gene/Prot IDs. This may take a few minutes...')
        fix_ids()
        print('Done.')