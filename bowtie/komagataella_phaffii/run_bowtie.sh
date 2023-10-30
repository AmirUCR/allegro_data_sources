#!/bin/bash

../bowtie/bowtie-build komagataella_phaffii_cds.fna index/komagataella_phaffii_cds
../bowtie/bowtie-build komagataella_phaffii_genomic.fna index/komagataella_phaffii_genomic

../bowtie/bowtie -v 0 -a index/komagataella_phaffii_cds komagataella_phaffii_a7_reads.fq output/a7/komagataella_phaffii_a7_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/komagataella_phaffii_cds komagataella_phaffii_a7_reads.fq output/a7/komagataella_phaffii_a7_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/komagataella_phaffii_cds komagataella_phaffii_a7_reads.fq output/a7/komagataella_phaffii_a7_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/komagataella_phaffii_cds komagataella_phaffii_e1_reads.fq output/e1/komagataella_phaffii_e1_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/komagataella_phaffii_cds komagataella_phaffii_e1_reads.fq output/e1/komagataella_phaffii_e1_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/komagataella_phaffii_cds komagataella_phaffii_e1_reads.fq output/e1/komagataella_phaffii_e1_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/komagataella_phaffii_genomic komagataella_phaffii_a7_reads.fq output/a7/komagataella_phaffii_a7_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/komagataella_phaffii_genomic komagataella_phaffii_a7_reads.fq output/a7/komagataella_phaffii_a7_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/komagataella_phaffii_genomic komagataella_phaffii_a7_reads.fq output/a7/komagataella_phaffii_a7_2mm_genomic_alignment.sam

../bowtie/bowtie -v 0 -a index/komagataella_phaffii_genomic komagataella_phaffii_e1_reads.fq output/e1/komagataella_phaffii_e1_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/komagataella_phaffii_genomic komagataella_phaffii_e1_reads.fq output/e1/komagataella_phaffii_e1_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/komagataella_phaffii_genomic komagataella_phaffii_e1_reads.fq output/e1/komagataella_phaffii_e1_2mm_genomic_alignment.sam