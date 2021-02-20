function BCFTOOLS_MERGE()
{
	cd $1
	ls *flt.vcf | cut -d'.' -f 1 | sort -u | while read id;do
	bcftools view -Oz -o ${id}.flt.vcf.gz ${id}.flt.vcf
	bcftools index ${id}.flt.vcf.gz
	done
	touch merge.sh
	echo  "bcftools merge \\"> merge.sh
	ls *flt.vcf.gz | while read id;do
	echo  "$id \\" >> merge.sh
	done
	echo "-o merged.bcftools.vcf" >> merge.sh
	bash merge.sh
}
