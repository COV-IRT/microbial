#make a list of the files
#use a bit of regex to cut out the prefixes
ls PE/*_1.fastq*|cut -f1 -d"_">PE_manifest;
#ls SE/*.fastq.gz>SE_manifest;



#make a manifest of the PE fastqs
for file in $(cat PE_manifest); do
    # build up a list of commands
    echo "./knead_workflow.sh -f "\"$file"_1.fastq.bz2\" -r "\"$file"_2.fastq.bz2\" -t 68 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt
done
#make a manifest of the SE fastqs
for file in $(cat SE_manifest); do
    # build up a list of commands
    echo "./knead_workflow.sh -f "$file".fastq -t 68 > /scratch/06176/jochum00/COVIRT19/stdout/${file}.log 2>&1" >> commands.txt