#make a list of the files
#cd /scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/Microbial_Pre-Processing/Step1_Trimmed_Filtered_Reads/Step1c_Trimmed_Filtered_Read_Data/;
ls *_1_reads_trim30_1.fq.gz>PE_manifest;
ls *_trim30.fq.gz>SE_manifest;
#use a bit of regex to cut out the prefixes
sed -i "s/_1_reads_trim30_1.fq.gz//g" manifest;

#make a manifest of the PE fastqs
for file in $(cat PE_manifest); do
    # build up a list of commands
    echo "./metaphlan3_workflow.sh -v -f $file\_1_reads_trim30_1.fq.gz -r $file\_1_reads_trim30_2.fq.gz -t 48 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
done
#make a manifest of the SE fastqs
for file in $(cat SE_manifest); do
    # build up a list of commands
    echo "./metaphlan3_workflow.sh -v -f $file\_trim30.fq.gz -t 48 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
done