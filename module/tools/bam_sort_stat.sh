
function BAM_SORT_STAT(){

	cd $1
	sorted_bam_dir=$1/sorted_bam
	SAMTOOLS_SORT $PWD
	SAMTOOLS_INDEX ${sorted_bam_dir}
	SAMTOOLS_FLAGSTAT ${sorted_bam_dir}
	SAMTOOLS_STAT ${sorted_bam_dir}
	SAMTOOLS_IDXSTAT ${sorted_bam_dir}

}


function BAM_STAT(){
	cd $1
	SAMTOOLS_INDEX $PWD
	SAMTOOLS_FLAGSTAT $PWD
	SAMTOOLS_STAT $PWD
	SAMTOOLS_IDXSTAT $PWD
}
