#Primeramente configuramos el directorio de trabajo en la carpeta que contiene nuestros datos
setwd("Users/MariaGarcia/Downloads/Bioinformatics Eli")
#Consiguientemente le damos nombre a nuestra sequencias que ya se encuentran agrupadas. En este caso opte por "Homo" ya que nuestras sequencias corresponden a una poblacion de 20 humanos.
Homos <- readDNAStringSet("sequences.fasta Midterm.fasta")
#Ahora aplicamos msa y asigmanos la nueva valiable (opte por la misma)
Homos <- msa(Homos)
#Con esta siguiente funcion podemos determinar la longitud de nuestra sequencia (resultado = 642)
nchar(Homos)
#Utilice esta funciones para observar mi alineamiento y asi poder completar el paso 2 de la tarea 
print(Homos, show="complete")
alFreq <- alphabetFrequency(Homos)
alFreq
# Homo_sapiens 6  es el sujeto que parece presentar la mayoria de mutaciones en los que los nucleotidos han sido reeplazados en varias localizaciones a traves del alineamiento.
#ran Blast for Homo 6 - Homo sapiens mutant hemoglobin beta chain (HBB) gene (100% adjession number)
Mutated <-readDNAStringSet("Homo6")

