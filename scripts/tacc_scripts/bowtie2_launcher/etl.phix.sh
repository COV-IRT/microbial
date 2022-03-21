###################PAIRED###########################
#create the necessary directories
mkdir bowtie2_out/Phix/unaligned
mkdir bowtie2_out/Phix/aligned
mkdir bowtie2_out/Phix/SAM
mkdir bowtie2_out/Phix/BAM

#make a manifest of the fastqs
for forward in bowtie2_out/GRCh38/unaligned/*.GRCh38.unaligned.trim30.fwd.fq.gz; do
    # check for the paired FASTQ file
    reverse=${forward%fwd.fq.gz}rev.fq.gz
    if [ -f $reverse ]; then
        # build up a list of commands
        #echo "./bowtie2-phix.sh -v -f $forward -r $reverse > /scratch/06176/jochum00/COVIRT19/stdout/${forward}.log 2>&1" >> commands.txt
	echo "./bowtie2-phix.sh -v -f $forward -r $reverse" >> commands.txt
    else
        echo "WARNING: File $forward appears to be missing its pair ($reverse)"
    fi
done



##############SINGLE####################
#create the necessary directories
mkdir bowtie2_out/Phix/unaligned
mkdir bowtie2_out/Phix/aligned
mkdir bowtie2_out/Phix/SAM
mkdir bowtie2_out/Phix/BAM

#make a manifest of the fastqs
for file in bowtie2_out/GRCh38/unaligned/*trim30.fq.gz; do
    # build up a list of commands
#    echo "./bowtie2-phix.sh -v -f $file > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
    echo "./bowtie2-phix.sh -v -f $file" >> commands.txt
done
