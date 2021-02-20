
################################################
################################################
#THIS SCRIPT DEFINE A CUFFLINKASSEMBLE FUNCTION#
################################################
################################################

function CUFFLINKS_ASSEMBLE(){

			cd $1
			source $INI_PATH/config/bin_config
			strand_str="fr-firststrand"

			ls *bam|cut -d"." -f 1 |sort -u |while read id;do

				cufflinks -g ${gtf} \
				          -b ${genome_fasta} \
				          --library-type ${strand_str} \
				          --max-multiread-fraction 0.25 \
				          --3-overhang-tolerance 2000 \
				          -o cuff_out_${id} \
				          -p ${thread} ${star_bam}
				          
				mv cuff_out_${id}/transcripts.gtf ${id}_transcripts.gtf

			done

}
