#!/bin/sh
conda run -n allegro python 1_make_proteome.py
conda run -n diamond ./2_create_diamond_db_for_all_species.sh
conda run -n diamond ./3_run_diamond_for_ref_species.sh
conda run -n allegro python 4_aggregate_orthogroups.py