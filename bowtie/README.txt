Please download this file and disable wrapping to preserve the spacings.

First run library_gen.ipynb files.
Then run run_bowtie.sh for each species inside its folder.
Then run find_offtargets.ipynb

Inside the folder of each species you will see:

Loose Files:
    species_a7_library.csv:
        Contains the guides that satisfy the A7 constraints for this species.
        This subset is sliced from the large library for 2434 species "a7_ns_to_drive.csv"
        There will be AT LEAST 7 guides in this file. There may be more.

    species_e1_library.csv:
        Contains the guides that satisfy the E1 constraints for this species.
        This subset is taken from the large library for 2434 species "e1_ns_to_drive.csv"
        There will be at least one guide for each ortholog in this species.
        For example, if this species has only 3 orthologs, you will have at least 3 guides.

    species_a7_reads.fq:
        The same guides from species_a7_library.csv. but in FASTQ form for bowtie usage.
        This file exists only for reference for "READ_X" ID look up.

    species_e1_reads.fq:
        Same as above, but for e1 guides.

    species_cds.fna:
        CDS file for this species.

    species_genomic.fna:
        Genome file for this species.

Directories:
    a7/
        a7_lib_offtargets.csv:
            Same as species_a7_library.csv, but augmented with 6 more columns:
            Column 0mm_cds_offtargets shows the number of exact-matching offtargets IN THE CDS file for this guide 
            Column 1mm_cds_offtargets shows the number of mismatch by exactly 1 (i.e. a SNP) offtargets IN THE CDS file for this guide.
            Column 2mm_cds_offtargets shows the number of mismatch by exactly 2 (2 base mutations) offtargets IN THE CDS file for this guide.

            Column 0mm_genome_offtargets shows the number of exact-matching offtargets IN THE GENOME file for this guide.
            Column 1mm_genome_offtargets shows the number of mismatch by exactly 1 (i.e. a SNP) offtargets IN THE GENOME file for this guide.
            Column 2mm_genome_offtargets shows the number of mismatch by exactly 2 (2 base mutations) offtargets IN THE GENOME file for this guide.

            To find where these offtargets occur, see the files below:

        cds_targets/
            guide_cds_targets.csv:
                Shows where this guide targets IN THE CDS file.
                NOTE:
                    The "aligned_seq" column shows the sequence of the guide when the strand is positive +
                    However, the sequence in this column will be the reverse complement of the guide when the strand is negative -
                    When the strand is negative, the "mismatch" column shows mutations in REVERSE. For example,

                    query_name  strand          reference_name          mapping_position         aligned_seq             mismatch            guide_w_pam                 guide           num_mutations
                    READ_1        -    lcl|NC_006071.1_cds_XP_503770.1_4087   1089          CCTTCTGTCATGAACGTCGTCAT  "5:G>C,17:G>T"   ATGACGACGTTCATGACAGAAGG   ATGACGACGTTCATGACAGA         2

                    "5:G>C,17:G>T" for this you find the 5th (0-based) base FROM the right (it will be a C) and swap it for G. Then you find the 17th (0-based) letter from
                    the right side of the aligned_seq (a T) and swap for G.
                    "5:G>C,17:G>T" means that the real sequence on the lcl|NC_006071.1_cds_XP_503770.1_4087 gene is actually CCTTCGGTCATGAACGTGGTCAT.
        
        genome_targets/
            guide_genome_targets.csv:
                Shows where this guide targets IN THE GENOME file. Above note applies.
