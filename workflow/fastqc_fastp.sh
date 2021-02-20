
function FASTQC_FASTP()
{

		cd $1

		if [ ! -d $2/raw_fastqc_report ];then
		    mkdir -p $2/raw_fastqc_report
		fi
		
		SOFTLINK_READS $PWD
		FASTQC_PRE $PWD
		mv *fastqc.html $2/raw_fastqc_report
		mv *fastqc.zip $2/raw_fastqc_report
		multiqc $2/raw_fastqc_report
		mv multiqc_* $2/raw_fastqc_report
		if [ ! -d ../clean_reads ];then
		    mkdir -p ../clean_reads
		fi 

		if skip_fastp;then
		    
		    fq_num=$(ls ../clean_reads/*clean.fq.gz | wc -l)
		    cd ../clean_reads 
		    tmp_clean_path=$PWD
		    if [ $fq_num -eq 0 ];then
		      echo "The clean fastq file suffix with clean.fq.gz does not exist in the $tmp_clean_path directory"
		      exit
		    else
		      echo "skipped FASTP,preparing start reads mapping "
		    fi
		else
		    echo "start FASTP running for raw reads in $raw_reads_dir"
		    cd $raw_reads_dir
		    FASTP_FILTER $PWD
		    if [ ! -d $2/fastp_report ];then
		        mkdir -p mkdir -p $2/fastp_report
		    fi 
		    mv $raw_reads_dir/*fastp.* $2/fastp_report
		fi 
		cd ../clean_reads

		echo "start fastqc  for clean reads in $PWD"
		FASTQC $PWD
		if [ ! -d $2/clean_fastqc_report ];then
		  mkdir -p $2/clean_fastqc_report
		fi 
		mv *fastqc.html $2/clean_fastqc_report
		mv *fastqc.zip $2/clean_fastqc_report
		multiqc $2/clean_fastqc_report
		mv multiqc_* $2/clean_fastqc_report
		echo "Completed running of reads quality control"

}
