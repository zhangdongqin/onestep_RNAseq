#!/bin/bash

################################################
################################################
## 2021_02_03  provided by zdq             #####   
## author:zhangdongqin2@126.com            #####
################################################
################################################

INI_PATH=$PWD
CONDA_PATH=""
. ./params.config
. ./bin/usage.sh
#conda_x=""`find ~/ -name conda | sort -u | head -n 1 | cut -d '/' -f -5`

#echo $conda_x
parseCommandLine() {
  # Special case that nothing was provided on the command line so print usage
  # - include this if it is desired to print usage by default
  if [ "$#" -eq 0 ]; then
    Usage
    exit 0
  fi
  # Indicate specification for single character options
  # - 1 colon after an option indicates that an argument is required
  # - 2 colons after an option indicates that an argument is optional, must use -o=argument syntax
  optstring="hi:o::v"
  # Indicate specification for long options
  # - 1 colon after an option indicates that an argument is required
  optstringLong="help,reads-path:,genome-fasta:,gtf:,hisat-index:,output-path::,coldata-file::,bam-dir::,salmon-index::,deseq-condition::,cpu::,reads-count::,conda_path::,skip-fastp::,skip-hisat::,skip-stringtie::,skip-gffcompare::,skip-bcftools::,skip-htseq::,skip-deseq::,deseq-output-prefix::,version"
  # - 2 colons after an option indicates that an argument is optional, must use --option=argument syntax
  # Parse the options using getopt command
  # - the -- is a separator between getopt options and parameters to be parsed
  # - output is simple space-delimited command line
  # - error message will be printed if unrecognized option or missing parameter but status will be 0
  # - if an optional argument is not specified, output will include empty string ''
  
  GETOPT_OUT=$(getopt --options $optstring --longoptions $optstringLong -- "$@")
  exitCode=$?
  if [ $exitCode -ne 0 ]; then
    echo ""
    Usage
    exit 1
  fi
  # The following constructs the command by concatenating arguments
  # - the $1, $2, etc. variables are set as if typed on the command line
  # - special cases like --option=value and missing optional arguments are generically handled
  #   as separate parameters so shift can be done below
  eval set -- "$GETOPT_OUT"
  # Loop over the options
  # - the error handling will catch cases were argument is missing
  # - shift over the known number of options/arguments
  while true; do
    #echo "Command line option is $opt"
    case "$1" in
      -h|--help) # -h or --help  Print usage
        Usage
        exit 0
        ;;
      -i|--reads-path) # -i inputFile or --input-file inputFile  Specify the input file
        # Input file must be specified so $2 can be used
        raw_reads_dir=$2

        if [ ! -d $raw_reads_dir ];then

              echo "the reads directory $raw_reads_dir not exist!"
              exit
            fi
        shift 2
        ;;

      -o|--output-path) # -o outputFile or --output-file outputFile  Specify the output file
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            output="./"
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            output=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      --genome-fasta)
        
        genome_fasta=$2

             if [ ! -f  $genome_fasta ];then              
                echo "the genome fasta file $genome_fasta does not exist"
                exit
             fi

        shift 2
        ;;
      --gtf)
        
        gtf=$2

             if [ ! -f  $gtf ];then              
                echo "the gtf file $gtf does not exist"
                exit
             fi

        shift 2
        ;;
      --hisat-index)
        
        index=$2

        if [ ! -f  ${index}.1.ht2 ];then              
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.2.ht2 ];then 
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.3.ht2 ];then 
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.4.ht2 ];then 
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.5.ht2 ];then 
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.6.ht2 ];then 
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.7.ht2 ];then 
           echo "the index file $index not exist"
           exit
        elif [ ! -f  ${index}.8.ht2 ];then 
           echo "the index file $index not exist"
           exit
        fi

        shift 2
        ;;
      --cpu)

        case "$2" in
          "")  # No output file so use default (check elsewhere)
            thread=8
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            thread=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;
      --bam-dir)
        
        bam_dir=$2
        shift 2
        ;;
      --salmon-index)        
        salmon_index=$2
        shift 2
        ;;
      --skip-fastp)
        
        case "$2" in
          "")  
            skip_fastp=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_fastp=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;
      --skip-hisat)
        
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            skip_hisat=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_hisat=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;
      --skip-stringtie)
        
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            skip_stringtie=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_stringtie=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;
      --skip-gffcompare)
        
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            skip_gffcompare=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_gffcompare=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      --skip-bcftools)
        
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            skip_bcftools=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_bcftools=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      --skip-htseq)
        
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            skip_htseq=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_htseq=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      --skip-deseq)
        
        case "$2" in
          "")  # No output file so use default (check elsewhere)
            skip_deseq=false
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            skip_deseq=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;  

      --conda_path)

        case "$2" in
          "")  # No output file so use default (check elsewhere)
            CONDA_PATH=`find ~/ -name conda | sort -u | head -n 1 | cut -d '/' -f -5`
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            CONDA_PATH=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;          

      --reads-count)

        case "$2" in
          "")  # No output file so use default (check elsewhere)
            reads_count=""
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            reads_count=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;      

      --coldata-file)

        case "$2" in
          "")  # No output file so use default (check elsewhere)
            col_data=""
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            col_data=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      --deseq-condition)

        case "$2" in
          "")  # No output file so use default (check elsewhere)
            comp_condition=""
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            comp_condition=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      --deseq-output-prefix)

        case "$2" in
          "")  # No output file so use default (check elsewhere)
            results_name="zdq_deseq2_results"
            shift 2  # Because output file is an empty string $2=''
            ;;
          *) # Output file has been specified so use it
            results_name=$2
            shift 2  # Because output file is $2
            ;;
        esac
        ;;

      -v|--version) # -v or --version  Print the version
        printVersion
        exit 0
        ;;
      --) # No more arguments
        shift
        break
        ;;
      *) # Unknown option - will never get here because getopt catches up front
        echo ""
        echo "Invalid option $1." >&2
        Usage
        exit 1
        ;;
    esac
  done

##help
##reads-path
##output-path 
##genome-fasta
##gtf
##hisat-index
##cpu
#coldata-file
#deseq-condition
#reads-count
#deseq-output-prefix
#conda_path
#skip-fastp
#skip-hisat
#skip-stringtie
#skip-gffcompare
#skip-bcftools
#skip-htseq
#skip-deseq 
#version ##/


  # Get a list of all command line options that do not correspond to dash options.
  # - These are "non-option" arguments.
  # - For example, one or more file or folder names that need to be processed.
  # - If multiple values, they will be delimited by spaces.
  # - Command line * will result in expansion to matching files and folders.
  shift $((OPTIND-1))
  additionalOpts=$*
}

program=$(basename $0)
version="1.1.0"
versionDate="2021-02-5"
CONDA_PATH=`find ~/ -name conda | sort -u | head -n 1 | cut -d '/' -f -5`

# Initialize variables
#inputFile=""
#outputFile=""
#additionalOpts=""
################################################
################################################
## print the required for program             ##
################################################
################################################
# Parse the command line options
# - pass all arguments to the function
parseCommandLine "$@"

# Print command line information

echo "   reads path:   $raw_reads_dir          "
echo "   genome ref:   $genome_fasta           "
echo "   genome gtf:   $gtf                    "
echo "   hisat2 idx:   ${index}.1.ht2          "
echo "   output dir:   $output                 "
echo "   conda path:   $CONDA_PATH             "
echo "   cpu core  :   $thread                 "
echo "   program -v:   $program--version:$version"

################################################
################################################
## include function files                     ##
################################################
################################################
. ./module/tools/bam_sort_stat.sh
. ./module/software/mapping/hisat.sh
. ./module/software/qc/fastqc.sh
. ./module/tools/softlinkFORunmatchedread.sh
. ./module/software/qc/fastp.sh
. ./module/software/statistic/samtools_sort.sh
. ./module/software/statistic/samtools_index.sh
. ./module/software/statistic/samtools_stat.sh
. ./module/software/statistic/samtools_flagstat.sh
. ./module/software/statistic/samtools_idxstat.sh
. ./workflow/fastqc_fastp.sh
. ./workflow/aligh_with_hisat.sh
. ./module/software/quantity/htseq.sh
. ./module/pipeline/count_matrix_build.sh
. ./module/software/assemble/stringtie.sh
. ./module/pipeline/stringtie_merge.sh
. ./workflow/stringtie_assemble_merge.sh
. ./module/software/variants/gffcompare.sh
. ./workflow/tmap2gtf_lnrna.sh
. ./workflow/bcftools_variants.sh
. ./module/tools/merge_vcf.sh
. ./module/software/variants/bcftools.sh
#cat ./config/function.config | while read id;do . $id;done
################################################
################################################
## FASTQC FOR RAW READS                       ##
################################################
################################################
if [  ! -d  ${output} ];then
      mkdir -p ${output}
fi

if $skip_fastp;then
    echo "SKIPPED FASTP RUNNING,PREPARING RUNNING MAPPING STEP...."
else 
    FASTQC_FASTP $raw_reads_dir $output
    cd $raw_reads_dir
    cd ../clean_reads
    clean_reads_dir=$PWD
fi

################################################
################################################
## HISAT        MAPPING                    #####
################################################
################################################

if $skip_hisat;then
    stringtie_work_dir=${bam_dir}
    echo "SKIPPED HSIAT2 RUNNING,PREPARING RUNNING NEXT STEP..."
else
    if [ ! -d ${output}/hisat_out ];then                  
          mkdir -p ${output}/hisat_out
          work_dir=${output}/hisat_out
    else 
          work_dir=${output}/hisat_out
    fi
    ALIGH_WITH_HISAT2 ${work_dir}
    if $salmon_quant;then
          stringtie_work_dir=${output}/hisat_out
    else
          stringtie_work_dir=${output}/hisat_out/sorted_bam
    fi
fi
################################################
################################################
##BCFTOOLS VARIANTS CALLING                #####
################################################
################################################
if [ ! -d $output/snp_indel ];then
     mkdir -p $output/snp_indel
fi

BCFTOOLS_VARIANTS ${output}/hisat_out/sorted_bam $output/snp_indel

if [ ! -f $output/snp_indel/merged.bcftools.vcf ];then
    echo "ERROR:THE MERGED VCF FILE $output/snp_indel/merged.bcftools.vcf DOES NOT EXIST"
    exit
fi


################################################
################################################
## STRINGTIE ASSEMBLING                    #####
################################################
################################################
if [  ! -d  ${output}/stringtie_out ];then                  
      mkdir -p ${output}/stringtie_out      
fi

STRINGTIE_ASSEMBLE_MERGE ${output}/hisat_out/sorted_bam ${output}/stringtie_out
if [  ! -d  ${output}/gffcompare ];then
      mkdir -p ${output}/gffcompare
fi
cd ${output}/gffcompare
GFF_COMPARE
mergedGTF=${output}/stringtie_out/stringtie_merged.gtf
TMAP2GTF ${output}/gffcompare ${output}/stringtie_out/gffcmp.stringtie_merged.gtf.tmap

################################################
################################################
## HTSEQ-COUNTING                             ##
################################################
################################################
echo "start htseq counting"
cd $output/hisat_out/sorted_bam
    
if [ ! -d $output/htseq_out ];then
    mkdir -p $output/htseq_out
fi 
HTSEQ_COUNT $PWD

################################################
################################################
## COUNT MATRIX BUILDING                      ##
################################################
################################################
echo "start expression count matrix building"
cd $output/htseq_out

if [ ! -d $output/deg_dir ];then
    mkdir -p $output/deg_dir
fi
COUNT_MATRIX_BUILD $PWD
################################################
################################################
####Differential gene analysis        ##########
################################################
################################################


