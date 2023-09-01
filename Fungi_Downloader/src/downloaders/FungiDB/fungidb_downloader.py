import os
import requests
import threading
import pandas as pd

from downloaders.FungiDB.add_gene_prot_names import fix_ids
from utils.path_generator import generate_dirs
from utils.name_processor import process_two_part_name, process_name


class FungiDB_Downloader:
    def __init__(self) -> None:
        fungidb = pd.read_csv('src/downloaders/FungiDB/fungidb.csv', sep='\t')
        fungidb = fungidb[fungidb['Is Reference Strain'] == 'yes'].reset_index(drop=True)

        fungidb['short_name'] = fungidb['Species'].apply(lambda x: process_two_part_name(x))

        manifest = pd.read_csv('src/downloaders/FungiDB/GenomeDataTypes_Summary.csv')
        manifest['short_name'] = manifest['Species'].apply(lambda x: process_two_part_name(x))

        self.shared = manifest[manifest['Species'].isin(fungidb['Species'])]
        self.shared = self.shared[self.shared['Is Reference Strain'] == 'yes']
        self.shared = self.shared[self.shared['Protein coding genes'] != 0]

        generate_dirs('data/FungiDB')

        self.output_file_name = 'data/FungiDB/fungidb_input_species.csv'


    def fetch_url(self, row, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, allocated_i):
        new_name = row['short_name']
        file_name = process_name(row['Species'])

        print(new_name)
        
        cds_file_name = f'data/FungiDB/cds/{file_name}_cds.fna'
        protein_file_name = f'data/FungiDB/proteomes/{file_name}.faa'
        genome_file_name = f'data/FungiDB/genomes/{file_name}_genomic.fna'

        protein_url = row['Protein Fasta Download Link']
        genome_url = row['Genome Fasta Download Link']
        cds_url = protein_url.replace('_AnnotatedProteins.fasta', '_AnnotatedCDSs.fasta')

        if not os.path.exists(cds_file_name) or not os.path.exists(protein_file_name) or not os.path.exists(genome_file_name):
            try:
                cds_r = requests.get(cds_url)
                if str(cds_r.content[0:15]) == "b'<!doctype html>'":
                    print('no cds')
                    return
                
                open(cds_file_name, 'wb').write(cds_r.content)

                prot_r = requests.get(protein_url)
                if str(prot_r.content[0:15]) == "b'<!doctype html>'":
                    print('no prot')
                    os.remove(cds_file_name)
                    return
                
                open(protein_file_name, 'wb').write(prot_r.content)

                genome_r = requests.get(genome_url)
                if str(genome_r.content[0:15]) == "b'<!doctype html>'":
                    print('no genome')
                    os.remove(cds_file_name)
                    os.remove(protein_file_name)
                    return
                
                open(genome_file_name, 'wb').write(genome_r.content)

            except Exception as e:
                print(f'Something went wrong while downloading {file_name}')
                print(e)

                if os.path.exists(f'{cds_file_name}.gz'):
                    os.remove(f'{cds_file_name}.gz')

                if os.path.exists(f'{genome_file_name}.gz'):
                    os.remove(f'{genome_file_name}.gz')

                if os.path.exists(f'{protein_file_name}.gz'):
                    os.remove(f'{protein_file_name}.gz')

                return

        names[allocated_i] = new_name
        genome_fasta_files[allocated_i] = f'{file_name}_genomic.fna'
        cds_fasta_files[allocated_i] = f'{file_name}_cds.fna'
        original_names[allocated_i] = file_name
        cds_urls[allocated_i] = f'wget {cds_url}'
        genome_urls[allocated_i] = f'wget {genome_url}'
        proteome_urls[allocated_i] = f'wget {protein_url}'


    def fetch_url_chunk(self, rows):
        threads = []

        names = [None] * len(rows)
        genome_fasta_files = [None] * len(rows)
        cds_fasta_files = [None] * len(rows)
        original_names = [None] * len(rows)
        cds_urls = [None] * len(rows)
        genome_urls = [None] * len(rows)
        proteome_urls = [None] * len(rows)

        i = 0
        for _, row in rows.iterrows():
            thread = threading.Thread(
                target=self.fetch_url,
                args=(row, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, i)
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

        chunk_size = 5
        chunks = [self.shared.iloc[i:i+chunk_size] for i in range(0, len(self.shared), chunk_size)]

        for url_chunk in chunks:
            names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls = self.fetch_url_chunk(url_chunk)
            all_names.extend(names)
            all_genome_fasta_files.extend(genome_fasta_files)
            all_cds_fasta_files.extend(cds_fasta_files)
            all_original_names.extend(original_names)
            all_cds_urls.extend(cds_urls)
            all_genome_urls.extend(genome_urls)
            all_proteome_urls.extend(proteome_urls)


        with open('data/FungiDB/names.txt', 'w') as f:
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