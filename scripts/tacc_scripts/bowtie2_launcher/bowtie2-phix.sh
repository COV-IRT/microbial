#!/bin/bash

# the following modules are required
# module load intel/17.0.4
# module load samtools/1.5
# module load bowtie/2.3.2

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THIS=$(basename $0)

HELP="
Usage: ${THIS} [OPTIONS]... -f <file>

Runs a RNAseq pipeline based on bowtie2 as part of the COVIRT project

required arguments:
  -f    FASTQ file containing forward reads (or unpaired reads)

optional arguments:
  -r    FASTQ file containing reverse reads for paired-end data
  -t    number of threads to use (default: 68)
  -h    show help message
  -v    verbose output
"

function usage() {
    echo "$HELP"
    exit 0
}


pairedEnd="false"
verbose="false"
threads="68"

while getopts ":hvf:r:p:" o; do
    case "${o}" in
    f) # forward read
        forwardReads=${OPTARG}
        ;;
    r) # reverse read
        reverseReads=${OPTARG}
        pairedEnd="true"
        ;;
    f) # number of threads to use
        threads=${OPTARG}
        ;;
    v) # verbose
        verbose="true"
        ;;
    h | *) # print help text
        usage
        ;;
    esac
done
shift $((OPTIND - 1))


# early error detection
# assert that bowtie2 is in $PATH
type bowtie2 > /dev/null 2>&1 || { echo "ERROR: bowtie2 not found" >&2 ; exit 1; }
# assert that samtools is in $PATH
type samtools > /dev/null 2>&1 || { echo "ERROR: samtools not found" >&2 ; exit 1; }
# assert that $forwardReads is a valid file
[ -f "$forwardReads" ] || { echo "ERROR: Forward reads file $forwardReads not found" >&2 ; exit 1; }
# assert that this is either a single ended run OR $reverseReads is a valid file
[[ "$pairedEnd" == "false" ]] || [ -f "$reverseReads" ] || { echo "ERROR: Reverse reads file $reverseReads not found" >&2 ; exit 1; }
# assert that $threads is a number
[[ "$threads" =~ ^[0-9]+$ ]] || { echo "ERROR: Number of threads specified was not a number ($threads)" >&2 ; exit 1; }

# print module and version information
if [[ "$verbose" == "true" ]]; then
    echo -n "LOG: " && module list 2>&1 | tail -n+2
    echo "LOG: bowtie2 version information:"
    echo ""
    bowtie2 --version
    echo ""
    echo "LOG: samtools version information:"
    echo ""
    samtools --version
    echo ""
fi


# start pipeline

fileBasename=$(basename $forwardReads .fq.gz)

if [[ "$pairedEnd" == "true" ]]; then
    echo "LOG: Starting paired-end pipeline at $(date -Iseconds)"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Forward Reads: $(realpath $forwardReads)"
        echo "LOG: Reverse Reads: $(realpath $reverseReads)"
    fi
    echo ""

    fileBasename=${fileBasename%_1*}

    echo "LOG: running bowtie2"
    echo ""

    bowtie2 \
        -x /scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Phix/phi_plus_SNPs \
        -1 $forwardReads \
        -2 $reverseReads \
        -p $threads \
        --very-sensitive \
        --un-conc-gz $PWD/bowtie2_out/Phix/unaligned/$fileBasename.fastq.gz \
        --al-conc-gz $PWD/bowtie2_out/Phix/aligned/$fileBasename.fastq.gz \
        -S $PWD/bowtie2_out/Phix/SAM/${fileBasename}.sam

    echo "LOG: job completed"

else # pairedEnd == "false"
    echo "LOG: Starting single-end pipeline at $(date -Iseconds)"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Reads: $(realpath $forwardReads)"
    fi
    echo ""

    echo "LOG: running bowtie2"
    echo ""

    bowtie2 \
        -x /scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Phix/phi_plus_SNPs \
        -U $forwardReads \
        -p $threads \
        --very-sensitive \
        --un-gz $PWD/bowtie2_out/Phix/unaligned/$fileBasename.fastq.gz \
        --al-gz $PWD/bowtie2_out/Phix/aligned/$fileBasename.fastq.gz \
        -S $PWD/bowtie2_out/Phix/SAM/${fileBasename}.SE.sam

fi
