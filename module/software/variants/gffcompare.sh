################################################
################################################
## THIS SCRIPT DEFINE A GFFCOMPARE FUNCTION#####
################################################
################################################

function GFF_COMPARE(){
		
		gffcompare -R \
		-r $gtf \
		$output/stringtie_out/stringtie_merged.gtf > log.txt
}


