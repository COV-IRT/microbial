##############SINGLE####################
#create the necessary directories
mkdir sbatches

#make a manifest of the fastqs
for file in *.cram; do
    # build up a list of commands
    echo "./cram2fastq.sh -v -f $file > /scratch/06176/jochum00/scratch/06176/jochum00/smith_dataset/WGS/stdout/${file}.log 2>&1" >> commands.txt
done
