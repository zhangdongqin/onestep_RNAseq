################################################
################################################
## THIS SCRIPT DEFINE RSEM MATRIX BUILD        #
################################################
################################################


function RSEM_MATRIX(){

			cd $1
			source $INI_PATH/config/bin_config

			touch rsem-generate-data-matrix.sh
			echo "$RSEM_MATRIX \\" > rsem-generate-data-matrix.sh
			
			ls -lh $output/rsem_out/*.genes.results | awk '{print $9}' | sort -u | while read id;
			do
			
			echo "$id \\" >> rsem-generate-data-matrix.sh
			
			done
			
			echo "> $output/genecounts.txt" >>　rsem-generate-data-matrix.sh

			if [ ! -f  $output/rsem_out/rsem-generate-data-matrix.sh ];then

				 echo "rsem-generate-data-matrix.sh file does not exist "
			
			else 
			
				 bash rsem-generate-data-matrix.sh
			
			fi

}

####$1 rsem calculation results directory needed for this function
