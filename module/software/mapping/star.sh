################################################
################################################
## THIS SCRIPT DEFINE A STAR   MAPPING FUNCTION#
################################################
################################################

function ALIG_STAR(){

			cd $1 
			source $INI_PATH/config/bin_config

			ls *gz|cut -d"." -f 1 |sort -u |while read id;
			do

			    if [ ! -f  ${id}.1.clean.fq.gz ];then               
			        echo "The ${id}.1.clean.fq.gz file  does not exist "
			      exit
			    elif [ ! -f  ${id}.2.clean.fq.gz ];then              
			        echo "The ${id}.2.clean.fq.gz file does not exist "
			        exit
			    fi
			    if [  -f  ${id}.1.clean.fq.gz ];then               
			         
			        read1=${id}.1.clean.fq.gz
			    elif [  -f  ${id}.2.clean.fq.gz ];then              
			         
			        read2=${id}.2.clean.fq.gz

			    fi
					
			STAR  --runThreadN ${thread}  \
					--twopassMode Basic --genomeDir ${star_index} \
					--readFilesIn ${read1} ${read2} \
					--readFilesCommand zcat \
					--outSAMtype BAM SortedByCoordinate \
					--chimSegmentMin 20 \
					--outFilterIntronMotifs RemoveNoncanonical \
					--outFilterMultimapNmax 20 \
					--alignIntronMin 20 \
					--alignIntronMax 1000000 \
					--alignMatesGapMax 1000000 \
					--outFilterType BySJout \
					--alignSJoverhangMin 8 \
					--alignSJDBoverhangMin 1 \
					--outFileNamePrefix ${id}

			done
}


