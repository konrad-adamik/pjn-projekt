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

# Analiza głównych składowych dla wcześniej utworzonej macierzy częstości
pca <- prcomp(dtmTfIdfBounds)

# Przygotowanie danych do wykresu
legendPCA <- paste(paste("d",1:20,sep = ""), rownames(dtmTfIdfBounds), sep = " => ")
xPCA <- pca$x[,1]
yPCA <- pca$x[,2]

# Eksport wykresu do pliku
plotFilePCA <- paste(
    outputDir,
    "pca.png",
    sep = "\\"
)
png(filename = plotFilePCA, width = 1920, height = 1080)

# Wykres w przestrzeni dwuwymiarowej
plot(
    xPCA,
    yPCA, 
    main = "Analiza g??wnych sk?adowych",
    xlab = "PC1",
    ylab = "PC2",
    col = "darkorchid4",
    #xlim = c(-0.06,0.06),
    #ylim = c(0,0.06),
    pch = 16
)
text(
    xPCA,
    yPCA,
    paste("d",1:20,sep = ""),
    col = "darkorchid4",
    pos = 2,
    cex = 0.8
)
legend(
    "bottom",
    legendPCA,
    cex = 1.2,
    text.col = "darkorchid4"
)

##-- Redukcja przy użyciu dekompozycji według wartości osobliwych --##

# Analiza ukrytych wymiarów semantycznych 
# (dekompozycja wg wartości osobliwych)

lsa <- lsa(tdmTfAllMatrix)

# Przygotowanie danych do wykresu
coordDocsLSA <- lsa$dk%*%diag(lsa$sk)
coordTermsLSA <- lsa$tk%*%diag(lsa$sk)
termsImportanceLSA <- diag(lsa$tk%*%diag(lsa$sk)%*%t(diag(lsa$sk))%*%t(lsa$tk))
importantTerms <- names(tail(sort(termsImportanceLSA),30))
coordImportantTermsLSA <- coordTermsLSA[importantTerms,]
ownTermsLSA <- c("pandemia")
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
#x2LSA <- coordOwnTermsLSA[,1]
#y2LSA <- coordOwnTermsLSA[,2]

# Eksport wykresu do pliku
plotFileLSA <- paste(
    outputDir,
    "lsa.png",
    sep = "\\"
)
png(filename = plotFileLSA, height = 1920, width = 1080)

# Wykres w przestrzeni dwuwymiarowej
plot(
    x1LSA,
    y1LSA, 
    main = "Analiza ukrytych wymiar?w semantycznych",
    xlab = "SD1",
    ylab = "SD2",
    col = "darkorchid4",
    #xlim = c(-25,5),
    #ylim = c(0,0.06),
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
    "topleft",
    legendLSA,
    cex = 1.2,
    text.col = "darkorchid4"
)


# Analiza skupień dokumentów (podpunkt E) ----------------------------------------------

##-- Metoda hierarchiczna --##
##-- Parametry metody --##
# 1. Macierz częstości
# a. Waga (weighting)
# b. Z zakresem uwzględnionych zmiennych (bounds)
# 2. Miara odległości (euclidean, jaccard, cosine)
# 3. Sposób wyznaczania odległości pomiędzy skupieniami 
#    (single, complete, ward.D2)

par(mai = c(1,2,1,1))
nDocumentsClust <- 20
legendClust <- paste(
    paste("d",1:20,sep = ""), 
    rownames(dtmTfAllMatrix), 
    sep = " => "
)
docNamesClust <- rownames(dtmTfAllMatrix)
rownames(dtmTfAllMatrix) <- paste("d",1:20,sep = "")
rownames(dtmTfBoundsMatrix) <- paste("d",1:20,sep = "")
patternClust <- c(3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,2,2,2,2,2)
coloursClust <- c("violet","orange","turquoise")
colorsHistClust <- c()
for (i in 1:nDocumentsClust){
    colorsHistClust[i] <- coloursClust[patternClust[i]]
}
names(colorsHistClust) <- paste("d",1:20, sep = "")

# Eksperyment 1
expDistMatrix1 <- dist(dtmTfAllMatrix, method = "euclidean")
hclust1 <- hclust(expDistMatrix1, method = "single")
plot(hclust1)
barplot(hclust1$height, names.arg = 19:1)

# Eksperyment 2
expDistMatrix2 <- dist(dtmTfAllMatrix, method = "cosine")
hclust2 <- hclust(expDistMatrix2, method = "ward.D2")
plot(hclust2)
barplot(hclust2$height, names.arg = 19:1)
expDendrogram2 <- as.dendrogram(hclust2)
coloredExpDendrogram2 <- color_branches(expDendrogram2, h = 1.2)
plot(coloredExpDendrogram2)
legend(
    "topright",
    legendClust,
    cex = 0.3,
)
coloredExpDendrogram2 <- color_branches(expDendrogram2, col = colorsHistClust[expDendrogram2 %>% labels])
plot(coloredExpDendrogram2)
nClusters <- 3
clusters2 <-cutree(hclust2, k = nClusters)
names(clusters2) <- docNamesClust
clustersMatrix2 <- matrix(0, nDocumentsClust, nClusters)
rownames(clustersMatrix2) <- docNamesClust
for (i in 1:nDocumentsClust) {
    clustersMatrix2[i,clusters2[i]] <- 1
}
corrplot(clustersMatrix2)

# Eksperyment 3
expDistMatrix3 <- dist(dtmTfBoundsMatrix, method = "jaccard")
hclust3 <- hclust(expDistMatrix3, method = "ward.D2")
plot(hclust3)
barplot(hclust3$height, names.arg = 19:1)
expDendrogram3 <- as.dendrogram(hclust3)
coloredExpDendrogram3 <- color_branches(expDendrogram3, h = 1.3)
plot(coloredExpDendrogram3)
legend(
    "topright",
    legendClust,
    cex = 0.3,
)
nClusters <- 3
clusters3 <-cutree(hclust3, k = nClusters)
names(clusters3) <- docNamesClust
clustersMatrix3 <- matrix(0, nDocumentsClust, nClusters)
rownames(clustersMatrix3) <- docNamesClust
for (i in 1:nDocumentsClust) {
    clustersMatrix3[i,clusters3[i]] <- 1
}
corrplot(clustersMatrix3)

# Porównanie wyników eksperymentów
Bk_plot(
    expDendrogram2,
    expDendrogram3,
    add_E = F,
    rejection_line_asymptotic = F,
    main = "Indeks Fawlkes'a Mallows'a",
    ylab = "Indeks Fawlkes'a Mallows'a"
)

randEx2Ex3 <- comPart(clusters2, clusters3)
randEx2Pattern <- comPart(clusters2, patternClust)
randEx3Pattern <- comPart(clusters3, patternClust)

# Eksperyment 4
nClusters4 <- 4
kmeans4 <- kmeans(dtmTfIdfBounds, centers = nClusters4)
clusters4 <- kmeans4$cluster
clustersMatrix4 <- matrix(0, nDocumentsClust, nClusters4)
rownames(clustersMatrix4) <- docNamesClust
for (i in 1:nDocumentsClust) {
    clustersMatrix4[i,clusters4[i]] <- 1
}
corrplot(clustersMatrix4)

randEx2Ex4 <- comPart(clusters2, clusters4)
randEx3Ex4 <- comPart(clusters3, clusters4)
randEx4Pattern <- comPart(clusters4, patternClust)


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

#-- Dla pierwszego dokumentu --#
#-- Waga tf jako miara ważności słów --#
keywordsTf1 <- head(sort(DTM_Tf_NoBounds_Matrix[1,],decreasing = TRUE))
keywordsTf1

# -- Waga tfidf jako miara ważności słów --#
keywordsTfIdf1 <- head(sort(dtmTfIdfBoundsMatrix[1,],decreasing = TRUE))
keywordsTfIdf1

#-- Prawdopodobieństwo w modelu LDA jako miara ważności słów --##
termsImportance1 <- c(experiment1$topics[1,]%*%experiment1$terms)
names(termsImportance1) <- colnames(experiment1$terms)
keywordsLda1 <- head(sort(termsImportance1,decreasing = TRUE))
keywordsLda1

#-- Chmura tagów --##
par(mai = c(0,0,0,0))
wordcloud(corpusPP[2], max.words = 200, colors = brewer.pal(8,"PuOr"))
