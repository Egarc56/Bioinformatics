# load libraries
library(msa)
library(Biostrings)
library(seqinr)

#First we configure the working directory in the folder that contains our data
setwd("/Users/MariaGarcia/Downloads/Bioinformatics Eli")
#Then I named the sequences that are already grouped. In this case I opted for "Homos" since the sequences correspond to a sample of 20 humans.
Homos <- readDNAStringSet("sequences.fasta Midterm.fasta")
#Now we apply msa and assign the new variable (opted for the same name)
Homos.aln <- msa(Homos) #rename new variables so that you don't erase the previous one
#With this following function we can determine the length of our sequence (result = 642)
nchar(Homos.aln)
#Use these functions to observe my alignment so I can complete step 2 of the Homework
print(Homos.aln, show="complete")
alFreq <- alphabetFrequency(Homos.aln)
alFreq
#Homo_sapiens 6 is the subject that appears to present mutations. Nucleotides have been replaced in various locations throughout the alignment.
#I ran Blast for Homo 6 
#Result- Homo sapiens mutant hemoglobin beta chain (HBB) gene
#Accession number-AY356351.1
#Here I am translating the Homo 6 sequence into protein
setwd("/Users/MariaGarcia/Downloads")
Mutant6 <- readDNAStringSet("Homosapien.fasta") # file does not exist in the Bioinformatics folder, so I can't access it

# another way to access that sample:
Mutant6 <- Homos$Homo_sapiens_6

MutantGene6 <- Biostrings::translate(Mutant6)
print(MutantGene6)
#Here I am printing the translated sequence into a fasta file using the write.fasta from the seqinr package which i saved to my Bioinformatics folder
write.fasta(sequences = MutantGene6, names = MutantGene6, file.out = "MutantHomoX.fasta") #make sure to provide file extension
# I used Blast and GenBank to search for the genetic diseases linked to this gene 
#After doing some research I found that this gene is associated with sickle beta 
# thalassaemia and sickle cell anemia. Yes, this individual has this genetic diseases

