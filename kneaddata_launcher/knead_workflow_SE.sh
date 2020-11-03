#!/bin/bash

# the following modules are required
module load tacc-singularity

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THIS=$(basename $0)

HELP="
Usage: ${THIS} [OPTIONS]... -f <file>

String of docker containers that collate and convert a cram
to paired end fastq files that go through a kneaddata pipeline

required arguments:

  -f    forward raw fastq.bz2 file
  -r    reverse raw fastq.bz2 file
optional arguments:
  -t    number of threads to use (default: 68)
  -m    the max-memory to be used (default:4g)
  -h    show help message

"
function usage() {
    echo "$HELP"
    exit 0
}



threads="68"
mem="4g"

while getopts ":hfr:t:m:" o; do
    case "${o}" in
    f) # forward fastq
        $forwardReads=${OPTARG}
        ;;
     f) # reverse fastq
        $reverseReads=${OPTARG}
        ;;
    t) # number of threads to use
        threads=${OPTARG}
        ;;
    m) # max memory
        mem=${OPTARG}
        ;;
    h | *) # print help text
        usage
        ;;
    esac
done
shift $((OPTIND - 1))

# assert that $threads is a number
[[ "$threads" =~ ^[0-9]+$ ]] || { echo "ERROR: Number of threads specified was not a number ($threads)" >&2 ; exit 1; }
fileBasename=$(basename $forwardReads _1.fastq.bz2)
# start pipeline
mkdir /scratch/06176/jochum00/COVIRT19/knead_out/$fileBasename;
#decompress to tmpdir for speed
lbzip2 -dkcv /scratch/06176/jochum00/COVIRT19/michael_files/raw/$forwardReads>$TMP/$fileBasename\_1.fastq;
lbzip2 -dkcv /scratch/06176/jochum00/COVIRT19/michael_files/raw/$reverseReads>$TMP/$fileBasename\_2.fastq;
kneaddata \
--input $TMP/$fileBasename\_1.fastq \
--input $TMP/$fileBasename\_2.fastq \
--run-fastqc-end \
--output-prefix $fileBasename \
--output /scratch/06176/jochum00/COVIRT19/knead_out/$fileBasename/ \
--threads $threads \
--processes $threads \
--trimmomatic /scratch/06176/jochum00/COVIRT19/reference_datasets/Trimmomatic-0.39 \
--reference-db /scratch/06176/jochum00/COVIRT19/reference_datasets/kneaddata/hg37 \
--log /scratch/06176/jochum00/COVIRT19/knead_out/$fileBasename/$fileBasename.knead.log \
-v \
--max-memory $mem

echo "LOG: kneaddata complete, cleaning tmpdir"
echo ""
#cleanup your tmpdir
rm -v  $TMP/$fileBasename\_1.fastq;
rm -v  $TMP/$fileBasename\_2.fastq;
#echo "rsync complete, exiting now"
echo "LOG: job completed"

exit 00