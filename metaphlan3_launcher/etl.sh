read_path="/scratch/06176/jochum00/COVIRT19/7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/Microbial_Pre-Processing/Step1_Trimmed_Filtered_Reads/Step1c_Trimmed_Filtered_Read_Data"


#make a manifest of the fastqs
for file in *trim25_1.fq.gz; do
    # build up a list of commands
    echo "./metaphlan3_workflow.sh -v -f $file > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> new_commands.txt
done
