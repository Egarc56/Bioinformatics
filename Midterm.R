#First we configure the working directory in the folder that contains our data
setwd("/Users/MariaGarcia/Downloads/Bioinformatics Eli")
#Then I named the sequences that are already grouped. In this case I opted for "Homos" since the sequences correspond to a sample of 20 humans.
Homos <- readDNAStringSet("sequences.fasta Midterm.fasta")
#Now we apply msa and assign the new variable (opted for the same name)
Homos <- msa(Homos)
#With this following function we can determine the length of our sequence (result = 642)
nchar(Homos)
#Use these functions to observe my alignment so I can complete step 2 of the Homework
print(Homos, show="complete")
alFreq <- alphabetFrequency(Homos)
alFreq
#Homo_sapiens 6 is the subject that appears to present mutations. Nucleotides have been replaced in various locations throughout the alignment.
#I ran Blast for Homo 6 
#Result- Homo sapiens mutant hemoglobin beta chain (HBB) gene
#Accession number-AY356351.1
#Here I am translating the Homo 6 sequence into protein
setwd("/Users/MariaGarcia/Downloads")
Mutant <- readDNAStringSet("sequence-7.fasta")
Mutant.seqinr <- Biostrings::translate(Mutant)
MutantGene <- Biostrings::translate(Mutant)
print(MutantGene)
#Here I am printing the translated sequence into a fasta file which i saved to my Bioinformatics folder
write.fasta(sequences = MutantGene, names = MutantGene, file.out = "Mutant")
# I ran Blastp to search the best match for my protein sequence 
#result: hemoglobin subunit beta isoform X1 [Mandrillus leucophaeus]
#Accession Number-XP_011830555.1
#After doing some research I found using GenBank that this gene is associated with sickle beta thalassaemia. Yes this individual has this genetic disease

