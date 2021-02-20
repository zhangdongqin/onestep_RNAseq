#!/bin/bash

bash hisat2_rnaseq.sh -t 16 -d /home/origene/yinshan.cui/rna-seq/reads/test_reads \
-g /home/data/genome_resequencing_DB/mus_musculus/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf \
-r /home/data/genome_resequencing_DB/mus_musculus/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
-i /home/origene/yinshan.cui/rna-seq/index/genome \
-o /home/origene/yinshan.cui/test/rna_output \
-c /home/origene/yinshan.cui/miniconda3
