#!/bin/bash

# the following modules are required
#!/bin/bash

# the following modules are required


DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THIS=$(basename $0)

HELP="
Usage: ${THIS} [OPTIONS]... -f <file>

Runs a metaphlan3 pipeline based on bowtie2 as part of the COVIRT project

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

while getopts ":hvf:r:t:" o; do
    case "${o}" in
    f) # forward read
        forwardReads=${OPTARG}
        ;;
    r) # reverse read
        reverseReads=${OPTARG}
        pairedEnd="true"
        ;;
    t) # number of threads to use
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
# assert that metaphlan is in $PATH
type metaphlan > /dev/null 2>&1 || { echo "ERROR: metaphlan not found" >&2 ; exit 1; }
# assert that $forwardReads is a valid file
[ -f "$forwardReads" ] || { echo "ERROR: Forward reads file $forwardReads not found" >&2 ; exit 1; }
# assert that this is either a single ended run OR $reverseReads is a valid file
[[ "$pairedEnd" == "false" ]] || [ -f "$reverseReads" ] || { echo "ERROR: Reverse reads file $reverseReads not found" >&2 ; exit 1; }
# assert that $threads is a number
[[ "$threads" =~ ^[0-9]+$ ]] || { echo "ERROR: Number of threads specified was not a number ($threads)" >&2 ; exit 1; }
# print module and version information
#if [[ "$verbose" == "true" ]]; then
#    echo -n "LOG: " && module list 2>&1 | tail -n+2
#    echo "LOG: metaphlan3 version information:"
#    echo ""
#    metaphlan --version
#    echo ""
#fi




# start pipeline
mkdir -p /tmp/metaphlan3/mpa_v30_CHOCOPhlAn_201901
mkdir -p $PWD/metaphlan3/bowtie2/
mkdir -p $PWD/metaphlan3/metaphlan3/profile/
mkdir -p $PWD/metaphlan3/sams/
mkdir -p $PWD/metaphlan3/bams/


# start pipeline

fileBasename=$(basename $forwardReads .fq.gz)

if [[ "$pairedEnd" == "true" ]]; then
    echo "LOG: Starting paired-end pipeline at $(date -Iseconds) on $fileBasename"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Forward Reads: $(realpath $forwardReads)"
        echo "LOG: Reverse Reads: $(realpath $reverseReads)"
    fi
    echo ""

    fileBasename=${fileBasename%_1*}

    echo "LOG: running metaphlan3 on $fileBasename"
    echo ""
# I removed this option to build the bowtie2db in the tmp dir because it takes a long time,
#but I'm not sure about the ability to 91 nodes to all read from the same db so I might put it back.
#        --bowtie2db /tmp/metaphlan3/mpa_v30_CHOCOPhlAn_201901 \
    metaphlan \
        $forwardReads,$reverseReads \
        --bowtie2db /scratch/06176/jochum00/COVIRT19/metaphlan3/mpa_v30_CHOCOPhlAn_201901 \
        --input_type fastq \
        --add_viruses \
        -s $PWD/metaphlan3/sams/$fileBasename.PE.sam \
        --bowtie2out $PWD/metaphlan3/bowtie2/$fileBasename.PE.bowtie2.bz2 \
        -o $PWD/metaphlan3/profiles/$fileBasename.PE.profile.tsv \
        --nproc $threads

    echo "LOG: running samtools sam to bam on $fileBasename"
    samtools view --threads $threads -bS $PWD/metaphlan3/sams/$fileBasename.PE.sam > $PWD/metaphlan3/bams/$fileBasename.PE.bam \
#    samtools view --threads $threads -f 4 -bS $PWD/metaphlan3/sams/$fileBasename.sam > $PWD/metaphlan3/bams/$fileBasename.unmapped.bam

    echo "LOG: $fileBasename job completed"

else # pairedEnd == "false"
    echo "LOG: Starting single-end pipeline at $(date -Iseconds) on $fileBasename"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Reads: $(realpath $forwardReads)"
    fi
    echo ""
    echo "LOG: running metaphlan3 on on $fileBasename"
    echo ""
    metaphlan \
            $forwardReads \
            --bowtie2db /tmp/metaphlan3/mpa_v30_CHOCOPhlAn_201901 \
            --input_type fastq \
            --add_viruses \
            -s $PWD/metaphlan3/sams/$fileBasename.SE.sam \
            --bowtie2out $PWD/metaphlan3/bowtie2/$fileBasename.SE.bowtie2.bz2 \
            -o $PWD/metaphlan3/profiles/$fileBasename.SE.profile.tsv \
            --nproc $threads

    echo "LOG: running samtools sam to bam on on $fileBasename"
    samtools view --threads $threads -bS $PWD/metaphlan3/sams/$fileBasename.SE.sam > $PWD/metaphlan3/bams/$fileBasename.SE.bam \
#   samtools view --threads $threads -f 4 -bS $PWD/metaphlan3/sams/$forwardReads.sam > $PWD/metaphlan3/bams/$forwardReads.unmapped.bam
    echo "LOG: $fileBasename job completed"
fi