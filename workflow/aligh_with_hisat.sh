
function ALIGH_WITH_HISAT2(){
		ALIGH_HISAT2 ${clean_reads_dir}
		if [  ! -f  ${output}/hisat_out/sorted_bam ];then
                	mkdir -p ${output}/hisat_out/sorted_bam
		fi
		if $salmon_quant;then 
			cd $1
			BAM_STAT $PWD
			if [ ! -d $output/bam_mapping_report ];then
		    	mkdir -p $output/bam_mapping_report
			fi
			mv *.stats $output/bam_mapping_report
			mv *.flagstat $output/bam_mapping_report
			mv *.idxstats $output/bam_mapping_report
		else 
			cd $1
			BAM_SORT_STAT $PWD
			if [ ! -d $output/sorted_bam_mapping_report ];then
		    	mkdir -p $output/sorted_bam_mapping_report
			fi
			mv *.stats $output/sorted_bam_mapping_report
			mv *.flagstat $output/sorted_bam_mapping_report
			mv *.idxstats $output/sorted_bam_mapping_report
		fi
}
