


################################################
################################################
#THIS SCRIPT DEFINE A STRINGTIE ASSEMBLE FUNCTION#
################################################
################################################

function STRINGTIE_ASSEMBLE(){

			cd $1
			ls *bam|cut -d"." -f 1 |sort -u |while read id;do
			  	if [ ! -f  ${id}.sorted.bam ];then               
			        echo "The ${id}.sorted.bam file  does not exist "
			      exit
			    fi

			    if [ -f  ${id}.sorted.bam ];then

			      hisat2_bam=${id}.sorted.bam
			    fi 
			stringtie \
				-p $thread \
				-G $gtf \
				-o $output/stringtie_out/${id}.transcripts.gtf \
				$hisat2_bam \
				-v \
				-A ${id}.genes_abundance.txt \
				-C ${id}.coverage.gtf \
				-b ${id}.ballgown
			done

}




