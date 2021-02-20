function COSEQ(){

		source $CONDA_PATH/etc/profile.d/conda.sh
		conda activate Rserver		
		Rscript $INI_PATH/scripts/coseq.r --count_file $1 --outdir $2
		conda deactivate 
}