

function Usage () {                 
    echo "Required parameters:"
    echo -e "\t--reads-path: /path/to/raw/reads/for/rnaseq/data"
    echo -e "\t--genome-fasta: genome fasta file for rnaseq analysis"
    echo -e "\t--gtf: genome annotation gtf file for rnaseq analysis"
    echo -e "\t--hisat-index: hisat2 genome index name for rnaseq data"
    echo "Optional parameters:"
    echo -e "\t--output-path: results output directoty for analysis,default is ./,you must specify output dir with the format of --output-path=/path/to/results"
    echo -e "\t--cpu:cpu cores specified for running program,default is 8core"
    echo -e "\t--bam-dir: bam file for stringtie analysis if you skipped hisat mapping"
    echo -e "\t--coldata-file: coldata file for deseq2 analysis"
    echo -e "\t--skip-fastp:skip reads quality control step for rnaseq,you can specify the parameter with '--skip-fastp=true'"
    echo -e "\t--skip-hisat:skip hisat2 mapping for the rnaseq analysis"
    echo -e "\t--skip-stringtie: skip transcripts assemble with stringtie for the pipeline"
    echo -e "\t--skip-gffcompare: skip transcripts compare"
    echo -e "\t--skip-bcftools: skip bcftools variants calling for the pipeline"
    echo -e "\t--skip-deseq: skip differential gene analysis for the pipeline"
    echo -e "\t--help: help information of the pipeline"
    echo -e "\t--salmon-index: transcripts salmon index path for salmon quantify"
    echo -e "\t--salmon-quant: optional quantify method for genes and transcripts quantify,you can use salmon with '--salmon-quant=true' in your analysis"
    exit -1
}
