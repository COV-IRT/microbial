#make a list of the files
ls /scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/Microbial_Pre-Processing/Step1_Trimmed_Filtered_Reads/Step1c_Trimmed_Filtered_Read_Data/*fq.gz>manifest;
#use a bit of regex to cut out the prefixes
sed -i "s/_filtered_1_reads_trim25_1.fq.gz//g" manifest;

#make a manifest of the fastqs
for file in $(cat manifest); do
    # build up a list of commands
    echo "./metaphlan3_workflow.sh -v -f $file -t 68 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
done
