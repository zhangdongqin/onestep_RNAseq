
################################################
################################################
## THIS SCRIPT DEFINE RSEM CALCULATION FUNCTION#
################################################
################################################

function RSEM_COUNT(){

			source $INI_PATH/config/bin_config

			ls *gz | cut -d'.' -f 1 | sort -u | while read id;do
			do
			    if [ ! -f  ${id}.1.clean.fq.gz ];then               
			        echo "The ${id}.1.clean.fq.gz file  does not exist "
			      	exit
			    fi
			    if [ ! -f  ${id}.2.clean.fq.gz ];then              
			        echo "The ${id}.2.clean.fq.gz file does not exist "
			        exit
			    fi
			    if [  -f  ${id}.1.clean.fq.gz ];then               
			         
			        read1=${id}.1.clean.fq.gz
			    fi
			    if [  -f  ${id}.2.clean.fq.gz ];then              
			         
			        read2=${id}.2.clean.fq.gz

			    fi 	

				$RSEM_CEXP \
				-p 8 \
				--bowtie2 \
				--bowtie2-path $1 \
				--paired-end \
				$read1 \
				$read2 \
				$2 \
				$output/rsem_out/${id}.bowtie2.rsem

			done

}




###$1 软件bowtie2的bin路径
###$2 rsem_bowtie2基因组索引