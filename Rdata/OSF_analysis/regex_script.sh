#########REGEX TO SEPERATE NUM AND FRACS###########
#nums
awk '{for(i=1;i<=9;i=i+1){printf "%s ", $i}{printf "%s", RS}}' bracken-species-combined-taxonomy.txt >tax
awk '{for(i=10;i<=NF;i=i+2){printf "%s ", $i}{printf "%s", RS}}' bracken-species-combined-taxonomy.txt >num
paste tax num>num_combined-bracken-outputs-with-lineages.tsv
#fracs
awk '{for(i=11;i<=NF;i=i+2){printf "%s ", $i}{printf "%s", RS}}' bracken-species-combined-taxonomy.txt >frac
paste tax frac>frac_combined-bracken-outputs-with-lineages.tsv
