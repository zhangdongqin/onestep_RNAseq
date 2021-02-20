################################################
################################################
## THIS SCRIPT DEFINE A FASTP FILTER FUNCTION  #
################################################
################################################

function FASTP_FILTER(){

		cd $1  
		ls *gz|cut -d"." -f 1 |sort -u |while read id;
		do
		  	if [ ! -f  ${id}.1.fq.gz ];then               
		         echo "The ${id}.1.fq.gz file  does not exist "
		       	 exit
		    fi

		    if [ ! -f  ${id}.2.fq.gz ];then             
		         echo "The ${id}.2.fq.gz file does not exist "
		         exit
		    fi    
		   
		    if [ -f  ${id}.1.fq.gz ];then               
		         
		     raw_read1=${id}.1.fq.gz

		    fi

		    if [ -f  ${id}.2.fq.gz ];then             
		         
		     raw_read2=${id}.2.fq.gz

		    fi    
		 	
			fastp \
			-i $raw_read1 \
			-o ../clean_reads/${id}.1.clean.fq.gz \
			-I $raw_read2 \
			-O ../clean_reads/${id}.2.clean.fq.gz    
			mv fastp.html ${id}.fastp.html
			mv fastp.json ${id}.fastp.json
		    
		    if [ ! -f  ../clean_reads/${id}.1.clean.fq.gz ];then              
		         echo "The fastp  for ../clean_reads/${id}.1.clean.fq.gzdoes not executed correctly "
		       	 exit
		    fi

		    if [ ! -f  ../clean_reads/${id}.2.clean.fq.gz ];then             
		         echo "The fastp  for ../clean_reads/${id}.2.clean.fq.gz does not executed correctly "
		         exit
		    fi    
		done

}


