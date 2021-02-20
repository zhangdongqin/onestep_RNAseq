
################################################
################################################
## THIS SCRIPT DEFINE A STRINGTIE MRGE FUNCTION#
################################################
################################################


function STRINGTIE_MERGE(){

		cd $1

		ls -lh $output/stringtie_out/*.transcripts.gtf| awk '{print $9}' > mergelist.txt
			if [ ! -s mergelist.txt ];then
				echo "the mergelist.txt is empty,please check"
				exit
			fi
		stringtie --merge \
		-p $thread \
		-G $gtf \
		-o $output/stringtie_out/stringtie_merged.gtf \
		mergelist.txt	

}


