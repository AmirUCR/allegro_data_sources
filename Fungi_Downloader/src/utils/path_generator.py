import os


def generate_dirs(base_dir):
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    cds_dir = os.path.join(base_dir, 'cds')
    genomes_dir = os.path.join(base_dir, 'genomes')
    prot_dir = os.path.join(base_dir, 'proteomes')

    if not os.path.exists(cds_dir):
        os.mkdir(cds_dir)

    if not os.path.exists(genomes_dir):
        os.mkdir(genomes_dir)
            
    if not os.path.exists(prot_dir):
        os.mkdir(prot_dir)
