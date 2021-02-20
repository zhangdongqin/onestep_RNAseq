################################################
################################################
## THIS SCRIPT DEFINE FASTQC  FUNCTION   #######
################################################
################################################

function FASTQC(){

			cd $1 

			ls *gz|cut -d"." -f 1 |sort -u |while read id;
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
				
			    fastqc -t $thread \
			    $read1 \
			    $read2

			done
}


function FASTQC_PRE(){

			cd $1 

			ls *gz|cut -d"." -f 1 |sort -u |while read id;
			do

			    if [ -f  ${id}.1.fq.gz ];then               
			         
			     raw_read1=${id}.1.fq.gz
			     
			    fi

			    if [ -f  ${id}.2.fq.gz ];then             
			         
			     raw_read2=${id}.2.fq.gz

			    fi   
				
			    fastqc -t $thread \
			    $raw_read1 \
			    $raw_read2

			done
}