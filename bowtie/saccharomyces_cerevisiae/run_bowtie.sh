#!/bin/bash

../bowtie/bowtie-build saccharomyces_cerevisiae_cds.fna index/saccharomyces_cerevisiae_cds
../bowtie/bowtie-build saccharomyces_cerevisiae_genomic.fna index/saccharomyces_cerevisiae_genomic

../bowtie/bowtie -v 0 -a index/saccharomyces_cerevisiae_cds saccharomyces_cerevisiae_a7_reads.fq output/a7/saccharomyces_cerevisiae_a7_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/saccharomyces_cerevisiae_cds saccharomyces_cerevisiae_a7_reads.fq output/a7/saccharomyces_cerevisiae_a7_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/saccharomyces_cerevisiae_cds saccharomyces_cerevisiae_a7_reads.fq output/a7/saccharomyces_cerevisiae_a7_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/saccharomyces_cerevisiae_cds saccharomyces_cerevisiae_e1_reads.fq output/e1/saccharomyces_cerevisiae_e1_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/saccharomyces_cerevisiae_cds saccharomyces_cerevisiae_e1_reads.fq output/e1/saccharomyces_cerevisiae_e1_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/saccharomyces_cerevisiae_cds saccharomyces_cerevisiae_e1_reads.fq output/e1/saccharomyces_cerevisiae_e1_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/saccharomyces_cerevisiae_genomic saccharomyces_cerevisiae_a7_reads.fq output/a7/saccharomyces_cerevisiae_a7_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/saccharomyces_cerevisiae_genomic saccharomyces_cerevisiae_a7_reads.fq output/a7/saccharomyces_cerevisiae_a7_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/saccharomyces_cerevisiae_genomic saccharomyces_cerevisiae_a7_reads.fq output/a7/saccharomyces_cerevisiae_a7_2mm_genomic_alignment.sam

../bowtie/bowtie -v 0 -a index/saccharomyces_cerevisiae_genomic saccharomyces_cerevisiae_e1_reads.fq output/e1/saccharomyces_cerevisiae_e1_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/saccharomyces_cerevisiae_genomic saccharomyces_cerevisiae_e1_reads.fq output/e1/saccharomyces_cerevisiae_e1_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/saccharomyces_cerevisiae_genomic saccharomyces_cerevisiae_e1_reads.fq output/e1/saccharomyces_cerevisiae_e1_2mm_genomic_alignment.sam