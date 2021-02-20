################################################
################################################
## THIS SCRIPT DEFINE A HTSEQ COUNT MATRIX BUILD
#				FUNCTION                      ##
################################################
################################################


function COUNT_MATRIX_BUILD(){

			cd $1

			ls *txt | cut -d'.' -f 1 |tr '\n' '\t' |awk '{print $0}'| sed 's/^/'\\t'&/g' > $output/deg_dir/counts.file
			paste *.txt | awk '{printf $1"\t";for(i=2;i<=NF;i+=2) printf $i"\t";printf $i"\n"}' > $output/deg_dir/reads_count.txt

			cd $output/deg_dir 
			cat reads_count.txt >> counts.file
			mv counts.file counts.txt

}


