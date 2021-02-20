################################################
################################################
## THIS SCRIPT DEFINE A STRINGTIE MRGE FUNCTION#
################################################
################################################


function STRINGTIE_ESTIMATE(){

		cd $1 

		ls *bam|cut -d"." -f 1 |sort -u |while read id;
		do
			if [ ! -f  ${id}.sorted.bam ];then               
		         echo "The ${id}.sorted.bam file  does not exist "
		     	 exit
		    fi
		    if [ ! -f  $output/stringtie_out/stringtie_merged.gtf ];then               
		         echo "The $output/stringtie_out/stringtie_merged.gtf file  does not exist "
		     	 exit
		    else 
		    	 stringtie_gtf=$output/stringtie_out/stringtie_merged.gtf
		    fi

		    if [ -f  ${id}.sorted.bam ];then

		      hisat2_bam=${id}.sorted.bam
		    fi 

			stringtie \
			-p $thread \
			$hisat2_bam \
			-G $stringtie_gtf \
			-e \
			-b $output/abundance_est/${id}
		done
}





