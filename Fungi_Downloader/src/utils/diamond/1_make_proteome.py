import os
import re
import sys
import yaml
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


if not os.path.exists('inputs'):
    os.makedirs('inputs/reference')
    os.makedirs('inputs/databases')

if not os.path.exists('inputs/reference'):
    os.mkdir('inputs/reference')

if not os.path.exists('inputs/databases'):
    os.mkdir('inputs/databases')

if not os.path.exists('outputs'):
    os.mkdir('outputs')

def map_gene_to_prot_id(
    cds_path: str, 
    gene_re=r'\[gene=(.*?)\]', 
    prot_id_re=r'\[protein_id=(.*?)\]'
    ) -> dict:

    map: dict[str, str] = dict()
    cds = list(SeqIO.parse(open(cds_path), 'fasta'))

    for record in cds:
        gene_match = re.search(gene_re, record.description)

        if gene_match:
            prot_id = re.search(prot_id_re, record.description)
            
            map[gene_match.group(1)] = prot_id.group(1) if prot_id else ''
    
    return map


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--config',
        type=argparse.FileType(mode='r'),
        default='make_proteome_config.yaml', 
        help='The config file to use. Must be placed in the root folder.',
    )

    args = parser.parse_args()
    if args.config:
        data = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__

        for key, value in data.items():
            arg_dict[key] = value

    return args
    

def main() -> int:
    args = parse_arguments()

    gene_names = args.gene_names
    
    gene_name_to_prot_id = map_gene_to_prot_id(args.cds_path)
    prot_id_to_gene_name: dict[str, str] = dict()

    for gene in gene_names:
        if gene in gene_name_to_prot_id:
            prot_id_to_gene_name[gene_name_to_prot_id[gene]] = gene

    proteome_records = list(SeqIO.parse(open(args.proteome_path), 'fasta'))

    seqs = list()
    for record in proteome_records:
        if record.id in prot_id_to_gene_name:

            sequence = SeqRecord(
                record.seq,
                id=record.id,
                description='[gene=' + prot_id_to_gene_name[record.id] + '] ' + record.description,
            )

            seqs.append(sequence)

    with open('inputs/reference/' + args.species + '.faa', 'w') as f:
        for seq in seqs:
            SeqIO.write(seq, f, 'fasta')

if __name__ == '__main__':
    sys.exit(main())