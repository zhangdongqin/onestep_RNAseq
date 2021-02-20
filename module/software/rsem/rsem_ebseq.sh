################################################
################################################
## THIS SCRIPT DEFINE RSEM EBSEQ FUNCTION      #
################################################
################################################


function RSEM_EBSEQ(){

			source $INI_PATH/config/bin_config
			source ~/.bashrc

			RSEM_RUN_EBSEQ \
				$output/genecounts.txt \
				$1 \
				rsem_EBSEQ.results

			RSEM_FDR_CONTROL \
				rsem_EBSEQ.results \
				0.05 rsem_EBSEQ.deg.txt 
}


### $1 sample matrix needed for this program

