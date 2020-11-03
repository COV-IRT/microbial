#make a list of the files
#cd /scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/Microbial_Pre-Processing/Step1_Trimmed_Filtered_Reads/Step1c_Trimmed_Filtered_Read_Data/;
ls PE/*_1.fastq*|cut -f1 -d"_">PE_manifest;
ls SE/*.fastq.gz>SE_manifest;
#use a bit of regex to cut out the prefixes
sed -i "s/_1_reads_trim30_1.fq.gz//g" PE_manifest;
sed -i "s/_trim30.fq.gz//g" SE_manifest;

#make a manifest of the PE fastqs
for file in $(cat PE_manifest); do
    # build up a list of commands
    echo "./knead_workflow.sh -v -f "$file"_1.fastq.bz2 -r "$file"_2.fastq.bz2 -t 68 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
done
#make a manifest of the SE fastqs
for file in $(cat SE_manifest); do
    # build up a list of commands
    echo "./knead_workflow.sh -v -f "$file".fastq -t 68 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt