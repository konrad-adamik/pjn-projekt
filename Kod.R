# Kod "wspólny" -----------------------------------------------------------


# Załadowanie bibliotek dla przetwarzania i stworzenia macierzy cz?sto?ci
library(tm)
library(hunspell)
library(stringr)

# Załadowanie biblioteki dla redukcji wymiarów przy użyciu dekompozycji według warto?ci osobliwych
library(lsa)

# Załadowanie bibliotek dla analizy skupień dokumentów
library(proxy)
library(dendextend)
library(corrplot)
library(flexclust)

# Załadowanie biblioteki dla analityki tematyk
library(topicmodels)

# Załadowanie biblioteki dla analizy słów kluczowych
library(wordcloud)

# Zmiana katalogu roboczego
workDir <- "D:\\Software\\UEK\\PJN_Projekt"  
setwd(workDir)

# Definicja ścieżek dostępu do katalogów funkcjonalnych
inputDir <- ".\\Dane"
outputDir <- ".\\Wyniki"

# Przetwarzanie wstępne (podpunkt A oraz B) ---------------------------------------------------


# Utworzenie korpusu dokumentów dla przetwarzania wstępnego
corpusPPDir <- paste(
    inputDir,
    "Dokumenty",
    sep = "\\"
)

corpusPP <- VCorpus(
    DirSource(
        corpusPPDir,
        "UTF-8",
        "*.txt"
    ),
    readerControl = list(
        language = "pl_PL"
    )
)

# Wstępne przetwarzanie
corpusPP <- tm_map(corpusPP, removeNumbers)
corpusPP <- tm_map(corpusPP, removePunctuation)
corpusPP <- tm_map(corpusPP, stripWhitespace)
corpusPP <- tm_map(corpusPP, content_transformer(tolower))

stoplistFilePP <- paste(
    inputDir,
    "stopwords_pl.txt",
    sep = "\\"
)
stoplistPP <- readLines(stoplistFilePP, encoding = "UTF-8")
corpusPP <- tm_map(corpusPP, removeWords, stoplistPP)
corpusPP <- tm_map(corpusPP, stripWhitespace)

# Usunięcie znaków em dash i 3/4
removeCharPP <- content_transformer(
    function(x,character) gsub(character, "",x)
)
corpusPP <- tm_map(corpusPP, removeCharPP, intToUtf8(8722))
corpusPP <- tm_map(corpusPP, removeCharPP, intToUtf8(190))

# Usunięcie rozszerzeń? z nazw dokumentów w korpusie
cutExtensionsPP <- function(document){
    meta(document,"id") <- gsub(
        pattern = "\\.txt$",
        replacement = "",
        meta(document,"id")
    )
    return(document)
}
corpusPP <- tm_map(corpusPP, cutExtensionsPP)

# Usunięcie podziału na akapity w dokumentach tekstowych
pasteParagraphsPP <- content_transformer(
    function(text, char) paste(text, collapse = char)
)
corpusPP <- tm_map(corpusPP, pasteParagraphsPP, " ")

# Lematyzacja
polishDictPP <- dictionary("pl_PL")

lemmatize <- function(text){
    simpleText <- str_trim(as.character(text))
    vectorizedText <- hunspell_parse(simpleText, dict = polishDictPP)
    
    lemmatizedText <- hunspell_stem(vectorizedText[[1]], dict = polishDictPP)
    for (i in 1:length(lemmatizedText)) {
        if(length(lemmatizedText[[i]]) == 0) lemmatizedText[i] <- vectorizedText[[1]][i]
        if(length(lemmatizedText[[i]]) > 1) lemmatizedText[i] <- lemmatizedText[[i]][1]
    }
    newText <- paste(lemmatizedText, collapse = " ")
    return(newText)
}
corpusPP <- tm_map(corpusPP, content_transformer(lemmatize))

# Eksport przetworzonego korpusu
preprocessedDir <- paste(
    outputDir,
    "Dokumenty_Przetworzone",
    sep = "\\"
)
dir.create(preprocessedDir, showWarnings = FALSE)
writeCorpus(corpusPP, path = preprocessedDir)


# Macierze częstości (podpunkt C) ------------------------------------------------------

# Funkcje pomocnicze do tworzenia macierzy częstości o określonych parametrach

createTDM <- function(weightingMode, lowerBound, upperBound){
  TDM <- TermDocumentMatrix(
    corpusPP,
    control = list(
      weighting = weightingMode,
      bounds = list(
        global = c(lowerBound, upperBound)
      )
    )
  )
  return(TDM)
}

createDTM <- function(weightingMode, lowerBound, upperBound){
  DTM <- DocumentTermMatrix(
    corpusPP,
    control = list(
      weighting = weightingMode,
      bounds = list(
        global = c(lowerBound, upperBound)
      )
    )
  )
  return(DTM)
}

# Funkcja pomocnicza do eksportu macierzy do pliku

exportMatrix <- function(matrix, fileName) {
  matrixFile <- paste(
    outputDir,
    fileName,
    sep = "\\"
  )
  write.table(
    matrix, 
    file = matrixFile, 
    sep = ";", 
    dec = ",", 
    col.names = NA
  )
}

# Tworzenie macierzy częstości

##-- TDM --##

TDM_Tf_NoBounds <- TermDocumentMatrix(corpusPP)
TDM_Tf_Bounds_3_14 <- createTDM(weightTf, 3, 14)
TDM_Tf_Bounds_4_20 <- createTDM(weightTf, 4, 20)
TDM_Tf_Bounds_2_18 <- createTDM(weightTf, 2, 18)

TDM_TfIdf_NoBounds <- TermDocumentMatrix(
  corpusPP,
  control = list(
    weighting = weightTfIdf
  )
)
TDM_TfIdf_Bounds_3_14 <- createTDM(weightTfIdf, 3, 14)
TDM_TfIdf_Bounds_4_20 <- createTDM(weightTfIdf, 4, 20)
TDM_TfIdf_Bounds_2_18 <- createTDM(weightTfIdf, 2, 18)

##-- DTM --##

DTM_Tf_NoBounds <- DocumentTermMatrix(corpusPP)
DTM_Tf_Bounds_3_14 <- createDTM(weightTf, 3, 14)
DTM_Tf_Bounds_4_20 <- createDTM(weightTf, 4, 20)
DTM_Tf_Bounds_2_18 <- createDTM(weightTf, 2, 18)

DTM_TfIdf_NoBounds <- DocumentTermMatrix(
  corpusPP,
  control = list(
    weighting = weightTfIdf
  )
)
DTM_TfIdf_Bounds_3_14 <- createDTM(weightTfIdf, 3, 14)
DTM_TfIdf_Bounds_4_20 <- createDTM(weightTfIdf, 4, 20)
DTM_TfIdf_Bounds_2_18 <- createDTM(weightTfIdf, 2, 18)


# Konwersje macierzy rzadkich do macierzy klasycznych

##-- TDM --##

TDM_Tf_NoBounds_Matrix <- as.matrix(TDM_Tf_NoBounds)
TDM_Tf_Bounds_3_14_Matrix <- as.matrix(TDM_Tf_Bounds_3_14)
TDM_Tf_Bounds_4_20_Matrix <- as.matrix(TDM_Tf_Bounds_4_20)
TDM_Tf_Bounds_2_18_Matrix <- as.matrix(TDM_Tf_Bounds_2_18)

TDM_TfIdf_NoBounds_Matrix <- as.matrix(TDM_TfIdf_NoBounds)
TDM_TfIdf_Bounds_3_14_Matrix <- as.matrix(TDM_TfIdf_Bounds_3_14)
TDM_TfIdf_Bounds_4_20_Matrix <- as.matrix(TDM_TfIdf_Bounds_4_20)
TDM_TfIdf_Bounds_2_18_Matrix <- as.matrix(TDM_TfIdf_Bounds_2_18)

##-- DTM --##

DTM_Tf_NoBounds_Matrix <- as.matrix(DTM_Tf_NoBounds)
DTM_Tf_Bounds_3_14_Matrix <- as.matrix(DTM_Tf_Bounds_3_14)
DTM_Tf_Bounds_4_20_Matrix <- as.matrix(DTM_Tf_Bounds_4_20)
DTM_Tf_Bounds_2_18_Matrix <- as.matrix(DTM_Tf_Bounds_2_18)

DTM_TfIdf_NoBounds_Matrix <- as.matrix(DTM_TfIdf_NoBounds)
DTM_TfIdf_Bounds_3_14_Matrix <- as.matrix(DTM_TfIdf_Bounds_3_14)
DTM_TfIdf_Bounds_4_20_Matrix <- as.matrix(DTM_TfIdf_Bounds_4_20)
DTM_TfIdf_Bounds_2_18_Matrix <- as.matrix(DTM_TfIdf_Bounds_2_18)


# Eksport macierzy częstości do pliku

##-- TDM --##

exportMatrix(TDM_Tf_NoBounds_Matrix, "TDM_Tf_NoBounds_Matrix.csv")
exportMatrix(TDM_Tf_Bounds_3_14_Matrix, "TDM_Tf_Bounds_3_14_Matrix.csv")
exportMatrix(TDM_Tf_Bounds_4_20_Matrix, "TDM_Tf_Bounds_4_20_Matrix.csv")
exportMatrix(TDM_Tf_Bounds_2_18_Matrix, "TDM_Tf_Bounds_2_18_Matrix.csv")
exportMatrix(TDM_TfIdf_NoBounds_Matrix, "TDM_TfIdf_NoBounds_Matrix.csv")
exportMatrix(TDM_TfIdf_Bounds_3_14_Matrix, "TDM_TfIdf_Bounds_3_14_Matrix.csv")
exportMatrix(TDM_TfIdf_Bounds_4_20_Matrix, "TDM_TfIdf_Bounds_4_20_Matrix.csv")
exportMatrix(TDM_TfIdf_Bounds_2_18_Matrix, "TDM_TfIdf_Bounds_2_18_Matrix.csv")

##-- DTM --##

exportMatrix(DTM_Tf_NoBounds_Matrix, "DTM_Tf_NoBounds_Matrix.csv")
exportMatrix(DTM_Tf_Bounds_3_14_Matrix, "DTM_Tf_Bounds_3_14_Matrix.csv")
exportMatrix(DTM_Tf_Bounds_4_20_Matrix, "DTM_Tf_Bounds_4_20_Matrix.csv")
exportMatrix(DTM_Tf_Bounds_2_18_Matrix, "DTM_Tf_Bounds_2_18_Matrix.csv")
exportMatrix(DTM_TfIdf_NoBounds_Matrix, "DTM_TfIdf_NoBounds_Matrix.csv")
exportMatrix(DTM_TfIdf_Bounds_3_14_Matrix, "DTM_TfIdf_Bounds_3_14_Matrix.csv")
exportMatrix(DTM_TfIdf_Bounds_4_20_Matrix, "DTM_TfIdf_Bounds_4_20_Matrix.csv")
exportMatrix(DTM_TfIdf_Bounds_2_18_Matrix, "DTM_TfIdf_Bounds_2_18_Matrix.csv")

# Redukcje wymiarów (podpunkt D) -------------------------------------------------------


##-- Redukcja przy użyciu analizy głównych składowych --##

presentPca <- function(DTM, fileName){
  pca <- prcomp(DTM)
  
  legendPCA <- paste(paste("d",1:20,sep = ""), rownames(DTM), sep = " => ")
  xPCA <- pca$x[,1]
  yPCA <- pca$x[,2]
  
  plotFilePCA <- paste(
    outputDir,
    fileName,
    sep = "\\"
  )
  png(filename = plotFilePCA, width = 1000, height = 1300)
  
  plot(
    xPCA,
    yPCA, 
    main = "Analiza głównych składowych",
    xlab = "PC1",
    ylab = "PC2",
    col = "darkorchid4",
    pch = 16,
  )
  text(
    xPCA,
    yPCA,
    paste("d",1:20,sep = ""),
    col = "darkorchid4",
    pos = 2,
    cex = 0.75
  )
  legend(
    "top",
    legendPCA,
    cex = 0.9,
    text.col = "darkorchid4"
  )
  
  dev.off()
}

presentPca(DTM = DTM_Tf_NoBounds, fileName = "PCA_DTM_Tf_NoBounds.png")
presentPca(DTM = DTM_Tf_Bounds_2_18, fileName = "PCA_DTM_Tf_Bounds_2_18.png")
presentPca(DTM = DTM_Tf_Bounds_3_14, fileName = "PCA_DTM_Tf_Bounds_3_14.png")
presentPca(DTM = DTM_Tf_Bounds_4_20, fileName = "PCA_DTM_Tf_Bounds_4_20.png")

presentPca(DTM = DTM_TfIdf_NoBounds, fileName = "PCA_DTM_TfIdf_NoBounds.png")
presentPca(DTM = DTM_TfIdf_Bounds_2_18, fileName = "PCA_DTM_TfIdf_Bounds_2_18.png")
presentPca(DTM = DTM_TfIdf_Bounds_3_14, fileName = "PCA_DTM_TfIdf_Bounds_3_14.png")
presentPca(DTM = DTM_TfIdf_Bounds_4_20, fileName = "PCA_DTM_TfIdf_Bounds_4_20.png")

##-- Redukcja przy użyciu dekompozycji według wartości osobliwych --##

# Analiza ukrytych wymiarów semantycznych 
# (dekompozycja wg wartości osobliwych)

presentLsa <- function(TDM, fileName, ownTerms) {

  lsa <- lsa(TDM)
  
  coordDocsLSA <- lsa$dk%*%diag(lsa$sk)
  coordTermsLSA <- lsa$tk%*%diag(lsa$sk)
  termsImportanceLSA <- diag(lsa$tk%*%diag(lsa$sk)%*%t(diag(lsa$sk))%*%t(lsa$tk))
  importantTerms <- names(tail(sort(termsImportanceLSA),30))
  coordImportantTermsLSA <- coordTermsLSA[importantTerms,]
  ownTermsLSA <- ownTerms
  coordOwnTermsLSA <- coordTermsLSA[ownTermsLSA,]
  
  legendLSA <- paste(
    paste("d",1:20,sep = ""), 
    rownames(coordDocsLSA), 
    sep = " => "
  )
  x1LSA <- coordDocsLSA[,1]
  y1LSA <- coordDocsLSA[,2]
  x2LSA <- coordImportantTermsLSA[,1]
  y2LSA <- coordImportantTermsLSA[,2]
  
  # Eksport wykresu do pliku
  plotFileLSA <- paste(
    outputDir,
    fileName,
    sep = "\\"
  )
  png(filename = plotFileLSA, height = 1400, width = 1800)
  
  # Wykres w przestrzeni dwuwymiarowej
  plot(
    x1LSA,
    y1LSA, 
    main = "Analiza ukrytych wymiarów semantycznych",
    xlab = "SD1",
    ylab = "SD2",
    col = "darkorchid4",
    pch = 16
  )
  text(
    x1LSA,
    y1LSA,
    paste("d",1:20,sep = ""),
    col = "darkorchid4",
    pos = 4,
    cex = 0.8
  )
  points(
    x2LSA,
    y2LSA,
    col = "magenta",
    pch = 17
  )
  text(
    x2LSA,
    y2LSA,
    col = "magenta",
    labels = rownames(coordImportantTermsLSA),
    pos = 3,
    cex = 0.8
  )
  legend(
    "left",
    legendLSA,
    cex = 1.0,
    text.col = "darkorchid4"
  )
  dev.off()
}

ownTerms <- c("zdrowie", "covid", "przedsiębiorstwo", "oprogramowanie", "zakażenie", "pacjent", "ekonomia", "literatura", "biznes", "zarządzać")

presentLsa(TDM = TDM_Tf_NoBounds, fileName = "LSA_TDM_Tf_NoBounds.png", ownTerms)
presentLsa(TDM = TDM_Tf_Bounds_2_18, fileName = "LSA_TDM_Tf_Bounds_2_18.png", ownTerms)
presentLsa(TDM = TDM_Tf_Bounds_3_14, fileName = "LSA_TDM_Tf_Bounds_3_14.png", ownTerms)
# Przy użyciu tej macierzy nie działało przekazanie ownTerms, niestety nie wiem dlaczego
presentLsa(TDM = TDM_Tf_Bounds_4_20, fileName = "LSA_TDM_Tf_Bounds_4_20.png", c())

presentLsa(TDM = TDM_TfIdf_NoBounds, fileName = "LSA_TDM_TfIdf_NoBounds.png", ownTerms)
presentLsa(TDM = TDM_TfIdf_Bounds_2_18, fileName = "LSA_TDM_TfIdf_Bounds_2_18.png", ownTerms)
presentLsa(TDM = TDM_TfIdf_Bounds_3_14, fileName = "LSA_TDM_TfIdf_Bounds_3_14.png", ownTerms)
# Przy użyciu tej macierzy nie działało przekazanie ownTerms, niestety nie wiem dlaczego
presentLsa(TDM = TDM_TfIdf_Bounds_4_20, fileName = "LSA_TDM_TfIdf_Bounds_4_20.png", c())


# Analiza skupień dokumentów (podpunkt E) ----------------------------------------------

##-- Metoda hierarchiczna --##
##-- Parametry metody --##
# 1. Macierz częstości
# a. Waga (weighting)
# b. Z zakresem uwzględnionych zmiennych (bounds)
# 2. Miara odległości (euclidean, jaccard, cosine)
# 3. Sposób wyznaczania odległości pomiędzy skupieniami 
#    (single, complete, ward.D2)

##-- Przygotwanie --##

## Tu wpisać porównywane macierze do przygotowania oraz indeksu FM

clustMatrix1 <- DTM_Tf_NoBounds_Matrix
clustMatrix2 <- DTM_TfIdf_Bounds_4_20_Matrix

## metody dla matrix i hclust, do eksperymentów używamy tych samych dla obu macierzy
matrixMethod = "cosine"
hclustMethod = "complete"

## liczba skupień

nClusters <- 12


par(mai = c(1,2,1,1))
nDocumentsClust <- 20
legendClust <- paste(
    paste("d",1:20,sep = ""), 
    rownames(clustMatrix1), 
    sep = " => "
)
##-- nazwy zawsze takie same
docNamesClust <- rownames(DTM_Tf_NoBounds_Matrix)
rownames(clustMatrix1) <- paste("d",1:20,sep = "")
rownames(clustMatrix2) <- paste("d",1:20,sep = "")
patternClust <- c(3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,2,2,2,2,2)
coloursClust <- c("violet","orange","turquoise")
colorsHistClust <- c()
for (i in 1:nDocumentsClust){
    colorsHistClust[i] <- coloursClust[patternClust[i]]
}
names(colorsHistClust) <- paste("d",1:20, sep = "")

# Dendogram dla clustMatrix1
expClustMatrix1 <- dist(clustMatrix1, method = matrixMethod)
hclust1 <- hclust(expClustMatrix1, method = hclustMethod)

plot(hclust1)
barplot(hclust1$height, names.arg = 19:1)

expClustDendrogram1 <- as.dendrogram(hclust1)

clusters1 <-cutree(hclust1, k = nClusters)
names(clusters1) <- docNamesClust

clustersMatrix1 <- matrix(0, nDocumentsClust, nClusters)
rownames(clustersMatrix1) <- docNamesClust
for (i in 1:nDocumentsClust) {
    clustersMatrix1[i,clusters1[i]] <- 1
}
corrplot(clustersMatrix1)

# Dendogram dla clustMatrix2
expClustMatrix2 <- dist(clustMatrix2, method = matrixMethod)
hclust2 <- hclust(expClustMatrix2, method = hclustMethod)

plot(hclust2)
barplot(hclust2$height, names.arg = 19:1)

expClustDendrogram2 <- as.dendrogram(hclust2)

clusters2 <-cutree(hclust2, k = nClusters)
names(clusters2) <- docNamesClust

clustersMatrix2 <- matrix(0, nDocumentsClust, nClusters)
rownames(clustersMatrix2) <- docNamesClust
for (i in 1:nDocumentsClust) {
    clustersMatrix2[i,clusters2[i]] <- 1
}
corrplot(clustersMatrix2)

# Porównanie wyników dla clustMatrix1 oraz clustMatrix2
Bk_plot(
    expClustDendrogram1,
    expClustDendrogram2,
    add_E = F,
    rejection_line_asymptotic = F,
    main = "Indeks Fawlkes'a Mallows'a",
    ylab = "Indeks Fawlkes'a Mallows'a"
)



# Analiza tematyk (podpunkt F) ---------------------------------------------------------

# Analiza ukrytej alokacji Dirichlet'a

coloursLDA <- c("violet", "orange", "turquoise", "darkseagreen","blue", "red", "magenta", "yellow", "brown", "cyan", "grey","chocolate","beige","chartreuse","coral","aquamarine","purple","salmon","sienna","tomato")

#Prezentacja tematów
presentTopics <- function(numberOfTopics, resultsLDA){
  for(i in 1:numberOfTopics){
    plotFileLDA <- paste(
      outputDir,
      paste(paste("Temat", i, sep=" "),".png", sep=""),
      sep = "\\"
    )
    png(filename = plotFileLDA, width = 1000, height = 500)
    topic <- head(sort(resultsLDA$terms[i,],decreasing = TRUE), 10)
    par(mar = c(3,10,3,1))
    barplot(
      rev(topic),
      horiz = TRUE,
      las = 1,
      main = paste("Temat", i, sep=" "),
      xlab = "Prawdopodobieństwo",
      col = coloursLDA[i]
    )
    dev.off()
  }
}

#Prezentacja dokumentów
presentDocuments <- function(resultsLDA){
  for(i in 1:20){
    plotFileLDA <- paste(
      outputDir,
      paste(paste("Dokument", i, sep=" "),".png", sep=""),
      sep = "\\"
    )
    png(filename = plotFileLDA, width = 1000, height = 500)
    document <- resultsLDA$topics[i,]
    barplot(
      document,
      horiz = TRUE,
      las = 1,
      main = rownames(resultsLDA$topics)[i],
      xlab = "Prawdopodobieństwo",
      col = coloursLDA
    )
    dev.off()
  }
}

#Eksperymenetowanie z LDA
LDAExperiment <- function(DTM_Tf, numberOfTopics){
  nTermsLDA <- ncol(DTM_Tf)
  nTopicsLDA <- numberOfTopics
  lda <- LDA(
    DTM_Tf,
    k = numberOfTopics,
    method = "Gibbs",
    control = list(
      burnin = 2000,
      thin = 100,
      iter = 3000
    )
  )
  resultsLDA <- posterior(lda)
  presentTopics(numberOfTopics, resultsLDA)
  presentDocuments(resultsLDA)
  return(resultsLDA)
}

#Eksperyment 1, DTM bez granic, 4 tematy
experiment1 <- LDAExperiment(DTM_Tf_NoBounds, 4)
# Udział słów w tematach
options(scipen = 5)
words1 <- c("społeczny")
round(experiment1$terms[,words1],4)

words2 <- c("zdrowie")
round(experiment1$terms[,words2],4)

words3 <- c("proces")
round(experiment1$terms[,words3],4)

words4 <- c("życie")
round(experiment1$terms[,words4],4)

words5 <- c("praca")
round(experiment1$terms[,words5],4)

words6 <- c("osoba")
round(experiment1$terms[,words6],4)

#Eksperyment 2, DTM z granicami 2-18, 20 tematów
experiment2 <- LDAExperiment(DTM_Tf_Bounds_2_18, 20)
# Udział słów w tematach
options(scipen = 5)
words1 <- c("przedsiębiorstwo")
round(experiment2$terms[,words1],20)

words2 <- c("behawioralny")
round(experiment2$terms[,words2],20)

words3 <- c("społeczny")
round(experiment2$terms[,words3],20)

words4 <- c("praca")
round(experiment2$terms[,words4],20)

words5 <- c("folklor")
round(experiment2$terms[,words5],20)

words6 <- c("wirus","sarscov","epidemia")
round(experiment2$terms[,words6],20)

words9 <- c("życie")
round(experiment2$terms[,words9],20)

words10 <- c("literatura")
round(experiment2$terms[,words10],20)

words11 <- c("covid")
round(experiment2$terms[,words11],20)

words12 <- c("dorastać")
round(experiment2$terms[,words12],20)

words13 <- c("projekt", "zarządzać")
round(experiment2$terms[,words13],20)

#Eksperyment 3, DTM z granicami 3-14, 2 tematy
experiment3 <- LDAExperiment(DTM_Tf_Bounds_3_14, 2)

words1 <- c("rynek")
round(experiment3$terms[,words1],2)


# Analiza słów kluczowych (podpunkt G) ----------------------------------------------------------

#funkcja pomocnicza do eksperymentów na macierzach DTM z wagami Tf i TfIdf
weightExperiment <- function(DTM_Matrix){
  for(i in 1:20){
    print(paste("Dokument",i,sep=" "))
    print(head(sort(DTM_Matrix[i,],decreasing = TRUE)))
  }
}



#-- Waga tf jako miara ważności słów --#

#eksperyment 1, macierz Tf bez granic
weightExperiment(DTM_Tf_NoBounds_Matrix)

#eksperyment 2, macierz Tf z granicami 2-18
weightExperiment(DTM_Tf_Bounds_2_18_Matrix)

#eksperyment 3, macierz Tf z granicami 3-14
weightExperiment(DTM_Tf_Bounds_3_14_Matrix)

#eksperyment 4, macierz Tf z granicami 4-20
weightExperiment(DTM_Tf_Bounds_4_20_Matrix)



# -- Waga tfidf jako miara ważności słów --#

#eksperyment 5, macierz TfIdf bez granic
weightExperiment(DTM_TfIdf_NoBounds_Matrix)

#eksperyment 6, macierz TfIdf z granicami 2-18
weightExperiment(DTM_TfIdf_Bounds_2_18_Matrix)

#eksperyment 7, macierz TfIdf z granicami 3-14
weightExperiment(DTM_TfIdf_Bounds_3_14_Matrix)

#eksperyment 8, macierz TfIdf z granicami 4-20
weightExperiment(DTM_TfIdf_Bounds_4_20_Matrix)



#-- Prawdopodobieństwo w modelu LDA jako miara ważności słów --##

#eksperyment 9
for(i in 1:20){
  termsImportance1 <- c(experiment1$topics[i,]%*%experiment1$terms)
  names(termsImportance1) <- colnames(experiment1$terms)
  print(paste("Dokument",i,sep=" "))
  print(head(sort(termsImportance1,decreasing = TRUE)))
}

#-- Chmura tagów --##
par(mai = c(0,0,0,0))
wordcloud(corpusPP[2], max.words = 200, colors = brewer.pal(8,"PuOr"))
