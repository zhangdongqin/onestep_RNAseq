

################################################
################################################
## THIS SCRIPT DEFINE BCFTOOLS SNPCALL FUNCTION#
################################################
################################################

function HTSEQ_COUNT(){

		source $conda_bash
		conda activate htseq

		cd $1
		ls *bam | sort -u | while read id;
		do
		sample=$(basename ${id} .sorted.bam)
		htseq-count -f bam -s no ${id} $gtf> $output/htseq_out/${sample}.txt
		done
		conda deactivate 
}


