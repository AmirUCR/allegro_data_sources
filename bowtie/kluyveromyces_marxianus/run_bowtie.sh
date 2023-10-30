#!/bin/bash

../bowtie/bowtie-build kluyveromyces_marxianus_cds.fna index/kluyveromyces_marxianus_cds
../bowtie/bowtie-build kluyveromyces_marxianus_genomic.fna index/kluyveromyces_marxianus_genomic

../bowtie/bowtie -v 0 -a index/kluyveromyces_marxianus_cds kluyveromyces_marxianus_a7_reads.fq output/a7/kluyveromyces_marxianus_a7_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/kluyveromyces_marxianus_cds kluyveromyces_marxianus_a7_reads.fq output/a7/kluyveromyces_marxianus_a7_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/kluyveromyces_marxianus_cds kluyveromyces_marxianus_a7_reads.fq output/a7/kluyveromyces_marxianus_a7_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/kluyveromyces_marxianus_cds kluyveromyces_marxianus_e1_reads.fq output/e1/kluyveromyces_marxianus_e1_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/kluyveromyces_marxianus_cds kluyveromyces_marxianus_e1_reads.fq output/e1/kluyveromyces_marxianus_e1_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/kluyveromyces_marxianus_cds kluyveromyces_marxianus_e1_reads.fq output/e1/kluyveromyces_marxianus_e1_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/kluyveromyces_marxianus_genomic kluyveromyces_marxianus_a7_reads.fq output/a7/kluyveromyces_marxianus_a7_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/kluyveromyces_marxianus_genomic kluyveromyces_marxianus_a7_reads.fq output/a7/kluyveromyces_marxianus_a7_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/kluyveromyces_marxianus_genomic kluyveromyces_marxianus_a7_reads.fq output/a7/kluyveromyces_marxianus_a7_2mm_genomic_alignment.sam

../bowtie/bowtie -v 0 -a index/kluyveromyces_marxianus_genomic kluyveromyces_marxianus_e1_reads.fq output/e1/kluyveromyces_marxianus_e1_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/kluyveromyces_marxianus_genomic kluyveromyces_marxianus_e1_reads.fq output/e1/kluyveromyces_marxianus_e1_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/kluyveromyces_marxianus_genomic kluyveromyces_marxianus_e1_reads.fq output/e1/kluyveromyces_marxianus_e1_2mm_genomic_alignment.sam