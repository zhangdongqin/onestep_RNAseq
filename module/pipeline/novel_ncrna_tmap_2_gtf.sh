
################################################
################################################
## THIS SCRIPT DEFINE A TMAP2GTF FUNCTION   ####
################################################
################################################


function TMAP2GTF() {

	        # filtering novel lncRNA based on cuffmerged trascripts
	        awk '$3 =="x"||$3=="u"||$3=="i"{print $0}' $1 > novel.gtf.tmap
	        #   excluding length smaller than 200 nt
	        awk '$10 >200{print}' novel.gtf.tmap > novel.longRNA.gtf.tmap
	        #   extract gtf
	        awk '{print $5}' novel.longRNA.gtf.tmap |perl $INI_PATH/bin/extract_gtf_by_name.pl ${mergedGTF} - >novel.longRNA.gtf
	        awk '{if($3=="exon"){print $0}}' novel.longRNA.gtf > novel.longRNA.format.gtf 
	        perl $INI_PATH/bin/get_exoncount.pl novel.longRNA.format.gtf  > novel.longRNA.exoncount.txt
	        # gtf2gff3
	        #check whether required
	        # get fasta from gtf
	        gffread novel.longRNA.gtf -g ${genome_fasta} -w novel.longRNA.fa -W			


}




