# the following modules are required
#!/bin/bash


# the following modules are required
module load samtools

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THIS=$(basename $0)

HELP="
Usage: ${THIS} [OPTIONS]... -f <file>

Runs a porechop pipeline based as part of the COVIRT project

required arguments:
  -i    FASTQ file containing forward reads (or unpaired reads)

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

while getopts ":hvi:r:t:" o; do
    case "${o}" in
    i) # forward read
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
# assert that porechop is in $PATH
type porechop > /dev/null 2>&1 || { echo "ERROR: porechop not found" >&2 ; exit 1; }

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
mkdir -p $PWD/porechopped/
mkdir -p $PWD/nanofilt/
mkdir -p $PWD/minimap2/
mkdir -p $PWD/sniffles/
mkdir -p $PWD/bandage/
mkdir -p $PWD/vg/

# start pipeline

fileBasename=$(basename $forwardReads .fastq.gz)

if [[ "$pairedEnd" == "true" ]]; then
    echo "LOG: Starting paired-end pipeline at $(date -Iseconds) on $fileBasename"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Forward Reads: $(realpath $forwardReads)"
        echo "LOG: Reverse Reads: $(realpath $reverseReads)"
    fi
    echo ""

    fileBasename=${fileBasename%*}

else # pairedEnd == "false"
    echo "LOG: Starting single-end pipeline at $(date -Iseconds) on $fileBasename"
    if [[ "$verbose" == "true" ]]; then
        echo "LOG: Reads: $(realpath $forwardReads)"
    fi
    echo ""
    echo "LOG: running porechop on on $fileBasename"
    echo ""
    porechop \
            --format fastq \
            -i $forwardReads \
            -o  $PWD/porechopped/$fileBasename.fastq \
            -v 1 \
            --discard_middle 

#    unpigz -v porechopped/$forwardReads;

    NanoFilt -l 500 --headcrop 10 < $PWD/porechopped/$fileBasename.fastq > $PWD/nanofilt/$fileBasename.fastq
    pigz -v  $PWD/porechopped/$fileBasename.fastq;
    pigz -v  $PWD/nanofilt/$fileBasename.fastq;

    minimap2 -t $threads -k 28 -K 5G -a -o  $PWD/minimap2/$fileBasename.sam reference/NC_045512.2.fa  $PWD/nanofilt/$forwardReads
    echo "LOG: running samtools sam to bam on on $fileBasename"

    samtools view --threads $threads -b $PWD/minimap2/$fileBasename.sam > $PWD/minimap2/$fileBasename.bam 
    samtools sort -@$threads --verbosity 1 -o $PWD/minimap2/$fileBasename.sorted.bam minimap2/$fileBasename.bam
    samtools calmd -bAr -@$threads $PWD/minimap2/$fileBasename.sorted.bam reference/NC_045512.2.fa > $PWD/minimap2/$fileBasename.sorted.md.bam
    samtools index -@ $threads -b $PWD/minimap2/$fileBasename.sorted.md.bam  $PWD/minimap2/$fileBasename.sorted.md.bai 
    sniffles -t $threads -m $PWD/minimap2/$fileBasename.sorted.md.bam -v $PWD/sniffles/$fileBasename.vcf
    echo "LOG: $fileBasename job completed"
fi