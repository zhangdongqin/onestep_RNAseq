function SALMON_QUANT(){

		cd $1 
		ls *clean.fq.gz | sort -u | while read id;
		do
			sample=${basename $id .1.clean.fq.gz}
		    if [ ! -f  ${id} ];then               
		        echo "The ${id} file  does not exist "
		      	exit
		    fi
		    if [ ! -f  ${sample}.2.clean.fq.gz ];then              
		        echo "The ${sample}.2.clean.fq.gz file does not exist "
		        exit
		    fi
		    if [  -f  ${id} ];then               			         
		        read1=${id}
		    fi
		    if [  -f  ${sample}.2.clean.fq.gz ];then              			         
		        read2=${sample}.2.clean.fq.gz
		    fi				
			if [  !-d  ${output}/salmon_quant ];then                  
		        mkdir -p ${output}/salmon_quant
		    fi			    	
		    salmon quant \
		        --geneMap $gtf \
		        --threads $thread \
		        --index $salmon_index \
		        -1 ${read1} \
		        -2 ${read2} \
		        -o $sample			
		done
}

