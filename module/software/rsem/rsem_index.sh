################################################
################################################
## THIS SCRIPT DEFINE RSEM INDEX FUNCTION      #
################################################
################################################

function RSEM_IDNEX(){

		source $INI_PATH/config/bin_config
		
		$RSEM_INDEX \
		--gtf $gtf \
		--bowtie2 \
		--bowtie2-path $1 \
		-p $thread \
		$genome_fasta \
		$2

}



###$1 需要bowtie2的bin路径
###$2 需要bowtie2的索引路径和名字



