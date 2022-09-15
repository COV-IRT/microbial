#make a list of the files
#cd /scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/Microbial_Pre-Processing/Step1_Trimmed_Filtered_Reads/Step1c_Trimmed_Filtered_Read_Data/;
ls *fastq.gz>manifest;
#use a bit of regex to cut out the prefixes
sed -i "s/fastq.gz//g" manifest;

#make a manifest of the fastqs
for file in $(cat manifest); do
    # build up a list of commands
    echo "./porechop_workflow.sh -v -f "$file".fastq.gz -t 64 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
done