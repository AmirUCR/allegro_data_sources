#!/bin/bash

../bowtie/bowtie-build yarrowia_lipolytica_cds.fna index/yarrowia_lipolytica_cds
../bowtie/bowtie-build yarrowia_lipolytica_genomic.fna index/yarrowia_lipolytica_genomic

../bowtie/bowtie -v 0 -a index/yarrowia_lipolytica_cds yarrowia_lipolytica_a7_reads.fq output/a7/yarrowia_lipolytica_a7_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/yarrowia_lipolytica_cds yarrowia_lipolytica_a7_reads.fq output/a7/yarrowia_lipolytica_a7_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/yarrowia_lipolytica_cds yarrowia_lipolytica_a7_reads.fq output/a7/yarrowia_lipolytica_a7_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/yarrowia_lipolytica_cds yarrowia_lipolytica_e1_reads.fq output/e1/yarrowia_lipolytica_e1_0mm_cds_alignment.sam
../bowtie/bowtie -v 1 -a index/yarrowia_lipolytica_cds yarrowia_lipolytica_e1_reads.fq output/e1/yarrowia_lipolytica_e1_1mm_cds_alignment.sam
../bowtie/bowtie -v 2 -a index/yarrowia_lipolytica_cds yarrowia_lipolytica_e1_reads.fq output/e1/yarrowia_lipolytica_e1_2mm_cds_alignment.sam

../bowtie/bowtie -v 0 -a index/yarrowia_lipolytica_genomic yarrowia_lipolytica_a7_reads.fq output/a7/yarrowia_lipolytica_a7_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/yarrowia_lipolytica_genomic yarrowia_lipolytica_a7_reads.fq output/a7/yarrowia_lipolytica_a7_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/yarrowia_lipolytica_genomic yarrowia_lipolytica_a7_reads.fq output/a7/yarrowia_lipolytica_a7_2mm_genomic_alignment.sam

../bowtie/bowtie -v 0 -a index/yarrowia_lipolytica_genomic yarrowia_lipolytica_e1_reads.fq output/e1/yarrowia_lipolytica_e1_0mm_genomic_alignment.sam
../bowtie/bowtie -v 1 -a index/yarrowia_lipolytica_genomic yarrowia_lipolytica_e1_reads.fq output/e1/yarrowia_lipolytica_e1_1mm_genomic_alignment.sam
../bowtie/bowtie -v 2 -a index/yarrowia_lipolytica_genomic yarrowia_lipolytica_e1_reads.fq output/e1/yarrowia_lipolytica_e1_2mm_genomic_alignment.sam