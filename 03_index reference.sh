bwa index flycatcher_rearranged.fa

java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/references/flycatcher_rearranged.fa O=/home/jmanthey/references/flycatcher_rearranged.dict

samtools faidx flycatcher_rearranged.fa
