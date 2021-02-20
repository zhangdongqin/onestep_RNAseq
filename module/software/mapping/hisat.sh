
################################################
################################################
## THIS SCRIPT DEFINE A HISAT2 MAPPING FUNCTION#
################################################
################################################

function ALIGH_HISAT2(){

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
				if [  ! -d  ${output}/hisat_out ];then                  
			        mkdir -p ${output}/hisat_out
			    fi
				if [  ! -d  ${output}/hisat_out/unmapped_reads ];then                  
			        mkdir -p ${output}/hisat_out/unmapped_reads
			    fi			    	
				hisat2 \
					-p $thread \
					-x $index \
					-1 $read1 \
					-2 $read2 \
					--un-conc-gz ${output}/hisat_out/unmapped_reads/${id}.unmapped.fastq.gz \
					--met-stderr \
					--new-summary \
					--dta \
					--no-mixed \
					--no-discordant \
					--summary-file ${output}/hisat_out/${id}.hisat2.summary.log \
					| samtools view -bS -F 4 -F 8 -F 256 - > ${output}/hisat_out/${id}.bam
		        if [ -f ${output}/hisat_out/unmapped_reads/${id}.unmapped.fastq.1.gz ]; then
		            mv ${id}.unmapped.fastq.1.gz ${id}.unmapped_1.fastq.gz
		        fi
		        if [ -f ${output}/hisat_out/unmapped_reads/${id}.unmapped.fastq.2.gz ]; then
		            mv ${id}.unmapped.fastq.2.gz ${id}.unmapped_2.fastq.gz
		        fi				
			done
}

