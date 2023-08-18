# Allegro Data Sources

This repository contains instructions and scripts that allow you to download genomes, CDS, and proteomes for more than 2,000 fungal species used in ALLEGRO. It also contains instructions for running orthology analysis to find the orthologs to the seven genes of interest in _Saccharomyces cerevisiae_.

## Prerequisites
1. First, download Miniconda [https://docs.conda.io/en/main/miniconda.html](https://docs.conda.io/en/main/miniconda.html)
2. Clone this repository by either clicking on the green Code on the top right and clicking "Download ZIP," or downloading [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and then `$ git clone https://github.com/AmirUCR/allegro_data_sources.git` in your desired directory.
3. These installation instructions are for Ubuntu 20.04.6 LTS. Create a conda environment and activate it:

```
conda create -n allegro python=3.10 -y
conda activate allegro
```

4. Install the required Python packages.

```
pip install pandas
pip install biopython
pip install PyYAML
pip install requests
```

## Gathering the Data
The order you should visit these files is as follows:

1. NCBI
2. FungiDB
3. EnsemblFungi
4. MycoCosm
5. copy_db_union.py
6. fourdbs_concat/fix.ipynb

For steps 1 through 4, visit here: https://drive.google.com/drive/folders/1FSkpgUBtfJ4NcyftYYQicKJLKNVGUR9d?usp=drive_link 
After downloading data from each database, run the appropriate `2_add_gene_prot_names.py` to homogenize protein/CDS ID conventions.

Then run `copy_db_union.py` to merge all databases and place the unique species in directory `fourdb_concat`. Run all cells in `dourdbs_concat/fix.ipynb` to remove non-unicode characters (e.g., punctuations) from file names. This will also remove species with six or seven missing orthologs (we find these species in the next step after running diamond. After diamond is finished, we find these species in `diamond/analysis_missing_orthologs.ipynb`).


## Generating Orthologs

Create a conda environment with Python version 2.7 and activate it. Install diamond:

```
conda create -n diamond python=2.7
conda activate diamond
conda install -c bioconda diamond
```

Remaining steps:

7. diamond/0_just_run_this.sh
8. ortholog_finder

At this point, we are ready to run DIAMOND to find orthologs to the _Saccharomyces cerevisiae_ genes specified in `diamond/make_proteome_config.yaml` Open this yaml file and update the `cds_path` and `proteome_path` to point to the absolute paths of `fourdbs_concat/cds/saccharomyces_cerevisiae_cds.fna` and `fourdbs_concat/proteomes/saccharomyces_cerevisiae.faa` Proceed with step 7 above `$ sh 0_just_run_this.sh` which will automatically run the other four scripts and produce `Orthogroups.tsv`

Copy `Orthogroups.tsv` and place it in the directory `ortholog_finder`. Switch back to the `allegro` conda environment and proceed with step 8 above `$ python find_orthogroup.py` This script generates a final species list and an output folder which we input to ALLEGRO. The cds files in the output of this script contain the appropriate orthologs for each species. Copy the final input species csv file, along with the generated CDS as input to ALLEGRO.