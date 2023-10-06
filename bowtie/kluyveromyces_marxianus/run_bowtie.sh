#!/bin/bash

../bowtie/bowtie -v 0 index/km_cds kluyveromyces_marxianus_a7_reads.fq output/km_a7_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 index/km_cds kluyveromyces_marxianus_a7_reads.fq output/km_a7_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 index/km_cds kluyveromyces_marxianus_a7_reads.fq output/km_a7_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 index/km_genomic kluyveromyces_marxianus_a7_reads.fq output/km_a7_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 index/km_genomic kluyveromyces_marxianus_a7_reads.fq output/km_a7_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 index/km_genomic kluyveromyces_marxianus_a7_reads.fq output/km_a7_2mm_genomic_alignment.sam

../bowtie/bowtie -v 0 index/km_cds kluyveromyces_marxianus_e1_reads.fq output/km_e1_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 index/km_cds kluyveromyces_marxianus_e1_reads.fq output/km_e1_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 index/km_cds kluyveromyces_marxianus_e1_reads.fq output/km_e1_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 index/km_genomic kluyveromyces_marxianus_e1_reads.fq output/km_e1_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 index/km_genomic kluyveromyces_marxianus_e1_reads.fq output/km_e1_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 index/km_genomic kluyveromyces_marxianus_e1_reads.fq output/km_e1_2mm_genomic_alignment.sam