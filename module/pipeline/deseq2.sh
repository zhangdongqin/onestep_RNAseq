
################################################
################################################
## THIS SCRIPT DEFINE BCFTOOLS SNPCALL FUNCTION#
################################################
################################################

function DESEQ2(){

		source $INI_PATH/config/bin_config
		source $CONDA_PATH/etc/profile.d/conda.sh
		conda activate Rserver

		Rscript $INI_PATH/bin/deseq.r --count_file $1 --coldata_file $2 --sample_con $3 --outprefix $4
		conda deactivate 
}


