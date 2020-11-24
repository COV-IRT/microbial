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

while getopts ":hf:r:t:m:" o; do
    case "${o}" in
    f) # forward fastq
        forwardReads=${OPTARG}
        ;;
    r) # reverse fastq
        reverseReads=${OPTARG}
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
TMP="/tmp"
# start pipeline
mkdir /scratch/06176/jochum00/COVIRT19/knead_out/$fileBasename;
mkdir /tmp/$fileBasename;
#decompress to tmpdir for speed
echo "LOG: decompressing $fileBasename to tmpdir"

    fileBasename=${forwardReads%_1*}

lbzip2 -dkcv $forwardReads >"/tmp/"$fileBasename/$fileBasename"_1.fastq";
lbzip2 -dkcv $reverseReads > "/tmp/"$fileBasename/$fileBasename"_2.fastq";
echo "LOG: appending $fileBasename headers"
reformat.sh t=68 overwrite=T in="/tmp/"$fileBasename/$fileBasename"_1.fastq" in2="/tmp/"$fileBasename/$fileBasename"_2.fastq" out="/tmp/"$fileBasename"/reformat."$fileBasename"_1.fastq" out2="/tmp/"$fileBasename"/reformat."$fileBasename"_2.fastq" trd addslash;

echo "LOG: running kneadddata on $fileBasename "
kneaddata \
--input "/tmp/"$fileBasename"/reformat."$fileBasename"_1.fastq" \
--input "/tmp/"$fileBasename"/reformat."$fileBasename"_2.fastq" \
--run-fastqc-end \
--output-prefix $fileBasename \
--output /tmp/$fileBasename/ \
--threads $threads \
--processes $threads \
--trimmomatic /scratch/06176/jochum00/COVIRT19/reference_datasets/Trimmomatic-0.39 \
--reference-db /scratch/06176/jochum00/COVIRT19/reference_datasets/kneaddata/hg37 \
--log /tmp/$fileBasename/$fileBasename.knead.log \
-v \
--max-memory $mem

echo "LOG: kneaddata complete, cleaning tmpdir"
echo ""
#cleanup your tmpdir
rm -v "/tmp/"$fileBasename/$fileBasename"_1.fastq";
rm -v "/tmp/"$fileBasename/$fileBasename"_2.fastq";
lbzip2 -v /tmp/$fileBasename/*fastq;
rsync --progress --recursive /tmp/$fileBasename/* /scratch/06176/jochum00/COVIRT19/knead_out/$fileBasename/
#echo "rsync complete, exiting now"
rm -rf /tmp/$fileBasename/ ;
echo "LOG: job completed"

exit 00