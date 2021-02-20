################################################
################################################
## THIS SCRIPT DEFINE BCFTOOLS SNPCALL FUNCTION#
################################################
################################################

function BCF_SNPCALL(){

cd $1 

ls *.bam|cut -d"." -f 1 |sort -u |while read id;
do
    if [ ! -f  ${id}.sorted.bam ];then               
       echo "The ${id}.sorted.bam file  does not exist "
       exit
    fi
    if [ -f  ${id}.sorted.bam ];then               
        
        hisat2_bam=${id}.sorted.bam

    fi

    bcftools mpileup -Ou --fasta-ref $genome_fasta $hisat2_bam \
    | bcftools call -Ou -mv  \
    | bcftools norm -Ou -f $genome_fasta \
    | bcftools view  > $output/snp_indel/${id}.flt.vcf

done
}
