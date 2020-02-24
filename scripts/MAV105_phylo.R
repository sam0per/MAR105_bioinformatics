# from https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/

install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")
library(ape)
library(phangorn)
library(seqinr)

if (!requireNamespace("BiocManager", quietly=TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("msa")
library(msa)

# browseVignettes("msa")

# mySequenceFile <- system.file("examples", "exampleDNA.fasta", package="msa")
mySequenceFile <- "./data/MAR105_sequence.fasta"
mySequences <- readDNAStringSet(mySequenceFile)
mySequences

## ----doAlignment-----------------------------------------------
myFirstAlignment <- msa(mySequences)
myFirstAlignment

## ----showWholeWidth--------------------------------------------
print(myFirstAlignment, show="complete")

## ----IntegratePDF2---------------------------------------------
msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

## ----VisualizePDF,results='asis'-------------------------------
msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis",
               showNames="none", showLogo="none", askForOverwrite=FALSE)

## ----OtherAlgorithms-------------------------------------------
myClustalWAlignment <- msa(mySequences, "ClustalW")
myClustalWAlignment
myClustalOmegaAlignment <- msa(mySequences, "ClustalOmega")
myClustalOmegaAlignment
myMuscleAlignment <- msa(mySequences, "Muscle")
myMuscleAlignment

## ----helpPrint,eval=FALSE--------------------------------------
#  help("print,MsaDNAMultipleAlignment-method")

## ----printExamples---------------------------------------------
print(myFirstAlignment)
print(myFirstAlignment, show="complete")
print(myFirstAlignment, showConsensus=FALSE, halfNrow=3)
print(myFirstAlignment, showNames=FALSE, show="complete")

dna_phy = msaConvert(myClustalWAlignment, "phangorn::phyDat")
mt <- modelTest(dna_phy)
print(mt)
(best_mod = mt$Model[order(mt$BIC, decreasing = FALSE)[1]])
dna_dist <- dist.ml(dna_phy, model="F81")
my_UPGMA <- upgma(dna_dist)
my_NJ  <- NJ(dna_dist)
plot(my_UPGMA, main="UPGMA")
plot(my_NJ, main = "Neighbor Joining")
