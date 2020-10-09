#!/bin/bash

# the following modules are required
module load tacc-singularity
#singularity pull docker://cmirzayi/metaphlanwithconda
#!/bin/bash

# the following modules are required
module load tacc-singularity


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
type metaphlan > /dev/null 2>&1 || { echo "ERROR: metaphlan not found" >&2 ; exit 1; }
# assert that $forwardReads is a valid file
[ -f "$forwardReads" ] || { echo "ERROR: Forward reads file $forwardReads not found" >&2 ; exit 1; }
# assert that this is either a single ended run OR $reverseReads is a valid file
[[ "$pairedEnd" == "false" ]] || [ -f "$reverseReads" ] || { echo "ERROR: Reverse reads file $reverseReads not found" >&2 ; exit 1; }
# assert that $threads is a number
[[ "$threads" =~ ^[0-9]+$ ]] || { echo "ERROR: Number of threads specified was not a number ($threads)" >&2 ; exit 1; }

# print module and version information
if [[ "$verbose" == "true" ]]; then
    echo -n "LOG: " && module list 2>&1 | tail -n+2
    echo "LOG: metaphlan3 version information:"
    echo ""
    metaphlan3 --version
    echo ""
fi




# start pipeline
mkdir -p /tmp/metaphlan3/mpa_v30_CHOCOPhlAn_201901 
mkdir -p $PWD/metaphlan3/bowtie2/
mkdir -p $PWD/metaphlan3/metaphlan3/profile/
mkdir -p $PWD/metaphlan3/sams/
mkdir -p $PWD/metaphlan3/bams/


if [[ "$pairedEnd" == "true" ]]; then
    echo "LOG: Starting paired-end pipeline at $(date -Iseconds)"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Forward Reads: $(realpath $forwardReads)"
        echo "LOG: Reverse Reads: $(realpath $reverseReads)"
    fi

    echo "LOG: running metaphlan3"
    echo ""
    metaphlan \
        $forwardReads\_filtered_1_reads_trim25_1.fq.gz,$reverseReads\_filtered_1_reads_trim25_1.fq.gz \
        --bowtie2db /tmp/metaphlan3/mpa_v30_CHOCOPhlAn_201901 \
        --input_type fastq \
        -s $PWD/metaphlan3/sams/$forwardReads.sam \
        --bowtie2out $PWD/metaphlan3/bowtie2/$forwardReads.bowtie2.bz2 \
        -o $PWD/metaphlan3/profiles/$forwardReads\_profile.tsv \
        --nproc $threads 

    echo "LOG: running samtools sam to bam"
    samtools view --threads $threads -bS $PWD/metaphlan3/sams/$forwardReads.sam > $PWD/metaphlan3/bams/$forwardReads.bam \
#    samtools view --threads $threads -f 4 -bS $PWD/metaphlan3/sams/$forwardReads.sam > $PWD/metaphlan3/bams/$forwardReads.unmapped.bam

    echo "LOG: job completed"

else # pairedEnd == "false"
    echo "LOG: Starting single-end pipeline at $(date -Iseconds)"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Reads: $(realpath $forwardReads)"
    fi
    echo ""
    echo "LOG: running metaphlan3"
    echo ""
    metaphlan \
            $forwardReads\_filtered_1_reads_trim25_1.fq.gz \
            --bowtie2db /tmp/metaphlan3/mpa_v30_CHOCOPhlAn_201901 \
            --input_type fastq \
            -s $PWD/metaphlan3/sams/$forwardReads.sam \
            --bowtie2out $PWD/metaphlan3/bowtie2/$forwardReads.bowtie2.bz2 \
            -o $PWD/metaphlan3/profiles/$forwardReads\_profile.tsv \
            --nproc $threads 
    
    echo "LOG: running samtools sam to bam"
    samtools view --threads $threads -bS $PWD/metaphlan3/sams/$forwardReads.sam > $PWD/metaphlan3/bams/$forwardReads.bam \
#    samtools view --threads $threads -f 4 -bS $PWD/metaphlan3/sams/$forwardReads.sam > $PWD/metaphlan3/bams/$forwardReads.unmapped.bam

fi
