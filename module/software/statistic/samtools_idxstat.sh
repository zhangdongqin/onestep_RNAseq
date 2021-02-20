function SAMTOOLS_IDXSTAT(){

		cd $1 ###${output}/hisat_out 
		ls *bam | sort -u | while read id;
		do
			sample=$(basename $id .bam)
		    if [ ! -f  ${id} ];then               
		       echo "The ${id} file  does not exist "
		       exit
		    fi
		    if [ -f  ${id} ];then               	        
		       bam=${id}
		    fi			
			samtools idxstats $bam > ${bam}.idxstats
		done
}