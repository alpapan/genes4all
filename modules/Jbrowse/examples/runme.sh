#!/bin/bash
samtools view -S -b toy.sam -o toy.bam
samtools index toy.bam
../../bin/prepare-refseqs.pl -fasta IC413239AaEcon74.fsa
../../bin/bam-to-json.pl -bam toy.bam -track test -key testing

