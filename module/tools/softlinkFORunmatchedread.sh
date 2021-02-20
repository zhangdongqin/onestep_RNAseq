
function SOFTLINK_READS(){

	cd $1
	fq_gz_num2=$(ls *_1.fq.gz | wc -l)
	fq_gz_num3=$(ls *_R1.fq.gz | wc -l)
	fastq_gz_num=$(ls *_1.fastq.gz | wc -l)

	if [ $fq_gz_num2 -ne 0 ];then
	    
	    ls *_1.fq.gz | while read id;do
	      sample=$(basename $id _1.fq.gz)
	      ln -s $id ${sample}.1.fq.gz
	      ls -s ${sample}_2.fq.gz ${sample}.2.fq.gz
	    done
	fi

	#######create soft link for the reads with unmatch with 'sample.1.fq.gz sample.2.fq.gz '
	if [ $fq_gz_num3 -ne 0 ];then
	    
	    ls *_R1.fq.gz | while read id;do
	      sample=$(basename $id _R1.fq.gz)
	      ln -s $id ${sample}.1.fq.gz
	      ls -s ${sample}_R2.fq.gz ${sample}.2.fq.gz
	      done
	fi

	#######create soft link for the reads with unmatch with 'sample.1.fq.gz sample.2.fq.gz '
	if [ $fastq_gz_num -ne 0 ];then
	    
	    ls *_1.fastq.gz | while read id;do
	      sample=$(basename $id _1.fastq.gz)
	      ln -s $id ${sample}.1.fq.gz
	      ls -s ${sample}_2.fastq.gz ${sample}.2.fq.gz
	     done
	fi
	
}
