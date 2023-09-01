import re
import os
import sys
import shutil
import zipfile
import requests
import threading
import pandas as pd

from utils.path_generator import generate_dirs
from utils.name_processor import process_two_part_name, process_name


class NCBI_Downloader:
    def __init__(self) -> None:
        self.output_file_name = 'data/NCBI/ncbi_input_species.csv'

        # -------------------------READ API KEY --------------------------
        self.api_key = str()
        with open('src/downloaders/NCBI/ncbi_key.txt', 'r') as f:
            api_key = f.readline().strip()

        if not api_key:
            print('Couldn\'t read API key in src/downloaders/NCBI/ncbi_key.txt')
            sys.exit(1)
        # -----------------------------------------------------------------

        generate_dirs('data/NCBI')

        if not os.path.exists('src/downloaders/NCBI/genome_result.txt'):
            print('ERROR: src/downloaders/NCBI/genome_result.txt does not exist for NCBI downloader.')
            sys.exit(1)

        ncbi_query_file = ''
        with open('src/downloaders/NCBI/genome_result.txt', 'r') as f:
            ncbi_query_file = f.read()

        organism_regex = r'^\d+\.\s(.*)\n'
        self.matches = re.findall(organism_regex, ncbi_query_file, flags=re.MULTILINE)


    def get_taxon_id(self, species_name, api_key) -> str:
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


    def get_accession(self,taxon_id, api_key) -> str:
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


    def fetch_url(self, name, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, allocated_i):
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
        
        new_name = process_two_part_name(name)
        file_name = process_name(name)

        print(new_name)

        cds_file_name = f'data/NCBI/cds/{file_name}_cds.fna'
        protein_file_name = f'data/NCBI/proteomes/{file_name}.faa'
        genome_file_name = f'data/NCBI/genomes/{file_name}_genomic.fna'

        tax_id = self.get_taxon_id(species_name=name, api_key=self.api_key)  # GET TAXON
        accession = self.get_accession(taxon_id=tax_id, api_key=self.api_key)  # GET ACCESSION

        # TEST FOR EXISTENCE
        if not os.path.exists(cds_file_name) or not os.path.exists(genome_file_name) or not os.path.exists(protein_file_name):
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
        cds_urls[allocated_i] =  f"curl -H \"Accept: application/zip\" \"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download?include_annotation_type=CDS_FASTA\" --output {file_name}_cds_download.zip"
        genome_urls[allocated_i] =  f"curl -H \"Accept: application/zip\" \"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA\" --output {file_name}_genome_download.zip"
        proteome_urls[allocated_i] =  f"curl -H \"Accept: application/zip\" \"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download?include_annotation_type=PROT_FASTA\" --output {file_name}_prot_download.zip"


    def fetch_url_chunk(self, chunk):
        threads = []

        names = [None] * len(chunk)
        genome_fasta_files = [None] * len(chunk)
        cds_fasta_files = [None] * len(chunk)
        original_names = [None] * len(chunk)
        cds_urls = [None] * len(chunk)
        genome_urls = [None] * len(chunk)
        proteome_urls = [None] * len(chunk)

        for i, name in enumerate(chunk):
            thread = threading.Thread(
                target=self.fetch_url,
                args=(name, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, i)
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
        
        chunk_size = 15
        chunks = [self.matches[i:i+chunk_size] for i in range(0, len(self.matches), chunk_size)]

        for url_chunk in chunks:
            names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls = self.fetch_url_chunk(url_chunk)
            all_names.extend(names)
            all_genome_fasta_files.extend(genome_fasta_files)
            all_cds_fasta_files.extend(cds_fasta_files)
            all_original_names.extend(original_names)
            all_cds_urls.extend(cds_urls)
            all_genome_urls.extend(genome_urls)
            all_proteome_urls.extend(proteome_urls)


        with open('data/NCBI/names.txt', 'w') as f:
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