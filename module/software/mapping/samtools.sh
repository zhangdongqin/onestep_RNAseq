   
################################################
################################################
## THIS SCRIPT DEFINE SAMTOOLS STAT FUNCTION####
################################################
################################################


function SAMTOOLS_STAT(){

		cd $1 
		ls *bam|cut -d"." -f 1 |sort -u |while read id;
		do
		    if [ ! -f  ${id}.hisat.bam ];then               
		       echo "The ${id}.hisat.bam file  does not exist "
		       exit
		    fi
		    if [ -f  ${id}.hisat.bam ];then               
		        
		        hisat2_bam=${id}.hisat.bam

		    fi
			samtools index ${hisat2_bam}
			samtools stats ${hisat2_bam} > ${id}.stats
			samtools flagstat ${hisat2_bam} > ${id}.flagstat
			samtools idxstats ${hisat2_bam} > ${id}.idxstats

		done
}















