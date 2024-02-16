#First we configure the working directory in the folder that contains our data
setwd("Users/MariaGarcia/Downloads/Bioinformatics Eli")
#Then I named the sequences that are already grouped. In this case I opted for "Homos" since the sequences correspond to a sample of 20 humans.
Homos <- readDNAStringSet("sequences.fasta Midterm.fasta")
#Now we apply msa and assign us the new varidable (opt for the same one)
Homos <- msa(Homos)
#With this following function we can determine the length of our sequence (result = 642)
nchar(Homos)
#Use these functions to observe my alignment so I can complete step 2 of the task
print(Homos, show="complete")
alFreq <- alphabetFrequency(Homos)
alFreq
# Homo_sapiens 6 is the subject that appears to present the majority of mutations in which nucleotides have been replaced in various locations throughout the alignment.
#ran Blast for Homo 6 - Homo sapiens mutant hemoglobin beta chain (HBB) gene
#Accession number-AY356351.1
#Here I am translating the Homo 6 sequence into protein
Mutant <-("seqdump.txt")
names(seqs) <-("AY356351.1 Homo sapiens mutant") 
Mutant <- msa(seqs)
print(Mutantx)
#Here I am printing the translated sequence into a fasta file which i saved to my Bioinformatics folder
write.fasta(names="Mutantx", sequences=Mutantx, file.out="Mutant.fasta")
# I ran Blastp to search the best match for my protein sequence 
#result: hemoglobin subunit beta isoform X1 [Mandrillus leucophaeus]
#Accession Number-XP_011830555.1
#After doing some research I found using GenBank that this gene is associated with sickle beta thalassaemia. Yes this individual has this genetic disease

