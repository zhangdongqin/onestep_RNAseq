function SAMTOOLS_SORT(){

		cd $1 ###${output}/hisat_out 
		ls *.bam|cut -d"." -f 1 |sort -u |while read id;
		do
		    if [ ! -f  ${id}.bam ];then               
		       echo "The ${id}.bam file  does not exist "
		       exit
		    fi
		    if [ -f  ${id}.bam ];then               
		        
		        hisat2_bam=${id}.bam
		    fi
		if [  ! -f  ${output}/hisat_out/sorted_bam ];then                  
	        	mkdir -p ${output}/hisat_out/sorted_bam
	    	fi			
			samtools sort  \
				-@ $thread \
				-o ${output}/hisat_out/sorted_bam/${id}.sorted.bam \
				-T $id \
				$hisat2_bam
		done
}
