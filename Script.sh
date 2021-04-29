#!/bin/bash
# Tools 
export PATH="/fast-data/Lanqing/Tools/Hg38_Mapping_Tools/bin":$PATH
printf "Pipeline Started at \t"
date +"%T"
######## Initialize Variable ########
ID=""
ref_name="hg38"
THREAD=4
LIB="NextSeq"
START_FROM_SCRATCH="n"
DIR=$PWD
INTERVAL="/fast-data/Lanqing/Tools/Hg38_Mapping_Tools/coveredRegion_final_use.hg38.bed"
######## Various functions ########
err(){
  echo "$1.exiting";
  exit 1; # any non-0 exit code signals an error
}

ckFile(){
  if [ ! -s "$1" ]; then
    err "$2 File '$1' not found or is empty";
  fi
}
is_in_path() {
  if [ ! -x "$(command -v "$1")" ]; then
    echo "$1 is not in path, existing"
    exit 1;
  fi
}

is_in_path "run-bwamem"
is_in_path "bwa"
is_in_path "samtools"
is_in_path "k8"
is_in_path "seqtk"
is_in_path "mosdepth"
is_in_path "fastqc"
#is_in_path "samblaster"
GATK=$(type -P "gatk")

usage() {
    if [ -n "$1" ]; then echo $1 >&2; fi
    cat <<EOF >&2
Usage: $(basename $0) [options] <ProjectTable>
Runs Exome pipeline using GATK gatk-4.1.4.1/

Options:
  -o <string>  	Output Directory  (default: $DIR)
  -r <ref>  	refference genome version either hg19 or hg38 (default: $ref_name)
  -i <string> 	Interval list from exome capture regions (default: $INTERVAL)
  -g <GATK> 	path to gatk (default: $GATK)
  -t <INT>	Number of Thread to use (default: $THREAD)
  -l <string>	Library name: LB in Bam (default: $LIB)
  -s <y/n>	Start everything from scrath ?(default: $START_FROM_SCRATCH)
  -h           	produce this message
  -d          	debug mode, do not delete intermediate files
EOF
    exit 1
}

while getopts ":hd:o:r:i:g:t:l:s:" opt; do
  case $opt in
	h)usage;;
	d)debug=1;;
	o)DIR=$OPTARG;;
	r)ref_name=$OPTARG ;;
	i)INTERVAL=$OPTARG;;
	g)GATK_JAR=$OPTARG ;;
	t)THREAD=$OPTARG;;
	l)LIB=$OPTARG;;
	s)START_FROM_SCRATCH=$OPTARG;;	
	\?)usage "Invalid option: -$OPTARG";;
	:)usage "Option -$OPTARG requires an argument";;
  esac
done

out_folder=$(dirname $out_prefix)
prefix="${File%%.*}"

i=${File##*/}
j=${File2##*/}
fastqDir=${File%/*}
name="$( cut -d '_' -f -3 <<< "$i" )"
sample_name=${out_prefix##*/}

if [ ! -d "$out_folder" ]
then
  mkdir -p $out_folder
fi

if [ ! -d "$out_folder/$sample_name/" ]
then
  mkdir -p $out_folder/$sample_name/
fi

cd $out_folder/$sample_name/

ProjectTable=$1
ckFile "$ProjectTable" "ProjectTable"


