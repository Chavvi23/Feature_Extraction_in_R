library(Peptides)
library(peptider)

fileConn<-file("log.txt","a")
args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  write(("At least one argument must be supplied for input file"), fileConn, append = TRUE)
} 

if(length(args)>1){
  write(c("More than one argument supplied for input file"), fileConn, append = TRUE)
}

if (!file.exists(args[1])) {
  write(c("File Not Exists"), fileConn, append = TRUE)
}
d = read.csv(args[1]) #accepting input from command prompt

write(("Input file argument taken"), fileConn, append = TRUE)

write(("Writing in output file started"), fileConn, append = TRUE)

i=1
for (sequence in d$Peptide.Sequence) {
  aIndex = aIndex(sequence)

  Boman = boman(sequence)

  InstaIndex = instaIndex(sequence)

  PPetide = ppeptide(sequence, libscheme = "NNK", N=10^8)

  Number_of_Neighbors = getNofNeighbors(sequence)

  Charge = charge(sequence)

  Hydrophobicity = hydrophobicity(sequence)

  Total_Neighbors = getNofNeighbors(sequence)

  Codons = codons(sequence,libscheme = "NNN")

  Hmoment_1 = hmoment(sequence, angle=100)
  Hmoment_2 = hmoment(sequence, angle=180)

  pI_Murray = pI(sequence,pKscale="Murray")
  pI_Dawson = pI(sequence,pKscale="Dawson")

  Length = lengthpep(sequence)

  Target = d$Target[i]

  Moelcular_weight = mw(sequence)

  MassChargeRatio = mz(sequence)

  MSWHIM = as.numeric(unlist(mswhimScores(seq=sequence)))
  names(MSWHIM) <- paste("MSWHIM_",c(1:length(MSWHIM)),sep=' ')

  stScales = as.numeric(unlist(stScales(seq=sequence)))
  names(stScales) <- paste("stScales_",c(1:length(stScales)),sep=' ')

  tScales = as.numeric(unlist(tScales(seq=sequence)))
  names(tScales) <- paste("tScales_",c(1:length(tScales)),sep=' ')

  Kidera = as.numeric(unlist(kideraFactors(seq=sequence)))
  names(Kidera) <- paste("Kidera_",c(1:length(Kidera)),sep=' ')

  AAComp = as.numeric(unlist(aaComp(seq=sequence)))
  names(AAComp) <- paste("AAComp_",c(1:length(AAComp)),sep=' ')

  BlosumIndices = as.numeric(unlist(blosumIndices(seq=sequence)))
  names(BlosumIndices) <- paste("BlosumIndices_",c(1:length(BlosumIndices)),sep=' ')

  cruciani_Properties = as.numeric(unlist(crucianiProperties(seq=sequence)))
  names(cruciani_Properties) <- paste("cruciani_Properties_",c(1:length(cruciani_Properties)),sep=' ')

  Fasgai_Vectors = as.numeric(unlist(fasgaiVectors(sequence)))
  names(Fasgai_Vectors) <- paste("fasgaiVectors_",c(1:length(Fasgai_Vectors)),sep=' ')

  if (i == 1){
    allFeatures = data.frame(sequence, Length, Moelcular_weight,Codons,MassChargeRatio,Total_Neighbors,aIndex,Boman,InstaIndex,PPetide,Number_of_Neighbors,Charge,Hydrophobicity,Hmoment_1,Hmoment_2, pI_Murray, pI_Dawson)
    final = c(allFeatures,MSWHIM,stScales,tScales,Kidera,AAComp,BlosumIndices,cruciani_Properties,Fasgai_Vectors,data.frame(Target))
    write.table(final, "output.csv", sep = ",", row.names=F, col.names = T, quote = F, append = T)
  }
  else
  {
    allFeatures = data.frame(sequence,Length,Moelcular_weight,Codons,MassChargeRatio,Total_Neighbors,aIndex,Boman,InstaIndex,PPetide,Number_of_Neighbors,Charge,Hydrophobicity,Hmoment_1,Hmoment_2, pI_Murray, pI_Dawson)
    final = c(allFeatures,MSWHIM,stScales,tScales,Kidera,AAComp,BlosumIndices,cruciani_Properties,Fasgai_Vectors,data.frame(Target))
    write.table(final, "output.csv", sep = ",", row.names=F, col.names = F, append = T)
  }
  i = i+1
}

write(c("Done writing in output file"), fileConn, append = TRUE)
