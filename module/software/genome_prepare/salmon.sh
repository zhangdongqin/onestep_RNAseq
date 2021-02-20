
function SALMON_INDEX(){

	    ################
		if [ ! -d $output/salmon_index ];then
	    	mkdir -p $output/salmon_index
		fi
		###############
	    gffread \
	    -w $output/salmon_index/transcripts.fa.tmp \
	    -g ${genome_fasta} \
	    $gtf
	    ###############
	    grep '^>' ${genome_fasta} | cut -d ' ' -f 1 > $output/salmon_index/decoys.txt
	    awk '{ if(/^>/){print $1}else{print $0}}' $output/salmon_index/transcripts.fa.tmp > $output/salmon_index/salmon_input_transcripts.fa
	    cat $output/salmon_index/salmon_input_transcripts.fa ${genome_fasta} > $output/salmon_index/gentrome.fa
	    ###############
	    salmon \
	    index \
	    --threads $thread \
	    -t $output/salmon_index/gentrome.fa \
	    -d $output/salmon_index/decoys.txt \
	    -i $output/salmon_index/salmon	
}
