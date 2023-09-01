import os
import re
import gzip
import shutil
import requests
import threading
import pandas as pd

from downloaders.EnsemblFungi.add_gene_prot_names import fix_ids
from utils.path_generator import generate_dirs
from utils.name_processor import process_two_part_name, process_name


class EnsemblFungi_Downloader:
    def __init__(self) -> None:
        self.df = pd.read_csv('src/downloaders/EnsemblFungi/ensemblFungi.csv')
        self.df = self.df.drop(columns=self.df.columns[0])

        generate_dirs('data/EnsemblFungi')

        self.output_file_name = 'data/EnsemblFungi/ensembl_input_species.csv'


    def fetch_url(self, row, names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls, allocated_i):
        new_name = process_two_part_name(row['Species'])
        file_name = process_name(row['Species'])
        
        print(new_name)

        cds_file_name = f'data/EnsemblFungi/cds/{file_name}_cds.fna'
        protein_file_name = f'data/EnsemblFungi/proteomes/{file_name}.faa'
        genome_file_name = f'data/EnsemblFungi/genomes/{file_name}_genomic.fna'

        cds_url = row['cds_url']
        genome_url = row['dna_url']
        protein_url = row['prot_url']

        # -- CDS --
        if not os.path.exists(cds_file_name) or not os.path.exists(genome_file_name) or not os.path.exists(protein_file_name):
            try:
                cds_r = requests.get(cds_url)
                if str(cds_r.content[0:15]) == "b'<!doctype html>'":
                    print(f'{file_name} no cds')
                    return

                m = re.findall(r"<a href=\"([^?/]+)\"", str(cds_r.content))
                index = [idx for idx, s in enumerate(m) if '.gz' in s][0]

                cds_url = cds_url + m[index]
                file = requests.get(cds_url)

                open(f'{file_name}_cds.gz', 'wb').write(file.content)
                with gzip.open(f'{file_name}_cds.gz', 'rb') as f_in:
                    with open(cds_file_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(f'{file_name}_cds.gz')
                
                # -- DNA --
                dna_r = requests.get(genome_url)
                if str(dna_r.content[0:15]) == "b'<!doctype html>'":
                    print('no genome')
                    os.remove(f'{file_name}_cds.gz')
                    return

                m = re.findall(r"<a href=\"([^?/]+dna.toplevel.fa.gz)\"", str(dna_r.content))
                index = [idx for idx, s in enumerate(m) if '.gz' in s][0]

                genome_url = genome_url + m[index]
                file = requests.get(genome_url)

                open(f'{file_name}_dna.gz', 'wb').write(file.content)
                with gzip.open(f'{file_name}_dna.gz', 'rb') as f_in:
                    with open(genome_file_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(f'{file_name}_dna.gz')
            
                # -- PROT --
                prot_r = requests.get(protein_url)
                if str(prot_r.content[0:15]) == "b'<!doctype html>'":
                    print(f'{file_name} no prot')
                    os.remove(f'{file_name}_cds.gz')
                    os.remove(f'{file_name}_dna.gz')
                    return
                
                m = re.findall(r"<a href=\"([^?/]+pep.all.fa.gz)\"", str(prot_r.content))
                index = [idx for idx, s in enumerate(m) if '.gz' in s][0]

                protein_url = protein_url + m[index]
                file = requests.get(protein_url)

                open(f'{file_name}_prot.gz', 'wb').write(file.content)
                with gzip.open(f'{file_name}_prot.gz', 'rb') as f_in:
                    with open(protein_file_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(f'{file_name}_prot.gz')

            except Exception as e:
                print(f'Something went wrong while downloading {file_name}')
                print(e)

                if os.path.exists(f'{cds_file_name}.gz'):
                    os.remove(f'{cds_file_name}.gz')

                if os.path.exists(f'{genome_file_name}.gz'):
                    os.remove(f'{genome_file_name}.gz')

                if os.path.exists(f'{protein_file_name}.gz'):
                    os.remove(f'{protein_file_name}.gz')

                if os.path.exists(f'{file_name}_dna.gz'):
                    os.remove(f'{file_name}_dna.gz')

                if os.path.exists(f'{file_name}_cds.gz'):
                    os.remove(f'{file_name}_cds.gz')

                if os.path.exists(f'{file_name}_prot.gz'):
                    os.remove(f'{file_name}_prot.gz')

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

        chunk_size = 10
        chunks = [self.df.iloc[i:i+chunk_size] for i in range(0, len(self.df), chunk_size)]

        for url_chunk in chunks:
            names, genome_fasta_files, cds_fasta_files, original_names, cds_urls, genome_urls, proteome_urls = self.fetch_url_chunk(url_chunk)
            all_names.extend(names)
            all_genome_fasta_files.extend(genome_fasta_files)
            all_cds_fasta_files.extend(cds_fasta_files)
            all_original_names.extend(original_names)
            all_cds_urls.extend(cds_urls)
            all_genome_urls.extend(genome_urls)
            all_proteome_urls.extend(proteome_urls)


        with open('data/EnsemblFungi/names.txt', 'w') as f:
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