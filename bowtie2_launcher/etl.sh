###################PAIRED###########################
#create the necessary directories
mkdir bowtie2_out
mkdir bowtie2_out/unaligned
mkdir bowtie2_out/aligned
mkdir bowtie2_out/SAM
mkdir bowtie2_out/BAM
mkdir sbatches

#make a manifest of the fastqs
for forward in *_1.fq.gz; do
    # check for the paired FASTQ file
    reverse=${forward%_1.fq.gz}_2.fq.gz
    if [ -f $reverse ]; then
        # build up a list of commands
        echo "./bowtie2-workflow.sh -v -f $forward -r $reverse > /scratch/06176/jochum00/COVIRT19/stdout/${forward}.log 2>&1" >> new_commands.txt
	#echo "./bowtie2-workflow.sh -v -f $forward -r $reverse" >> commands.txt
    else
        echo "WARNING: File $forward appears to be missing its pair ($reverse)"
    fi
done



##############SINGLE####################
#create the necessary directories
mkdir bowtie2_out
mkdir bowtie2_out/unaligned
mkdir bowtie2_out/aligned
mkdir bowtie2_out/SAM
mkdir bowtie2_out/BAM
mkdir sbatches

#make a manifest of the fastqs
for file in *trim30.fq.gz; do
    # build up a list of commands
    echo "./bowtie2-workflow.sh -v -f $file > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> new_commands.txt
    #echo "./bowtie2-workflow.sh -v -f $file" >> commands.txt
done
