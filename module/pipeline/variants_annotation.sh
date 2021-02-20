################################################
################################################
## THIS SCRIPT DEFINE ANNOVAR ANNOTATION    ####
################################################
################################################


function VARIANTS_ANNTATION(){

      source $INI_PATH/config/bin_config
      cd $1
      ls *raw.vcf | cut -d'.' -f 1 | sort -u | while read id;do
      bcftools view -Oz -o $output/gatk_vcf/${id}.raw.vcf.gz ${id}.raw.vcf
      bcftools index $output/gatk_vcf/${id}.raw.vcf.gz
      done

      touch merge.sh

      echo  "bcftools merge \\" > merge.sh
      ls *raw.vcf.gz | while read id;do
      echo  "$id \\" >> merge.sh
      done
      echo "-o merge.raw.all.vcf" >> merge.sh
      if [ -f  merge.sh ];then               
          bash merge.sh
      else
          echo "$output/gatk_vcf/merge.sh does not exist!"
          exit
      fi
      $CONVERT2ANNOVAR -format vcf4 merged.raw.all.vcf > variants.avinput
      $TABLE_ANNOVAR variants.avinput $humandb -buildver hg38 -out variants_annotation \
      -remove -protocol refGene,cytoBand,esp6500siv2_all \
      -operation g,r,f -nastring . -csvout      

}




