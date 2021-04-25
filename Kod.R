# Kod "wspólny" -----------------------------------------------------------


# Za³adowanie bibliotek dla przetwarzania i stworzenia macierzy czêstoœci
library(tm)
library(hunspell)
library(stringr)

# Za³adowanie biblioteki dla redukcji wymiarów przy u¿yciu dekompozycji wed³ug wartoœci osobliwych
library(lsa)

# Za³adowanie bibliotek dla analizy skupieñ dokumentów
library(proxy)
library(dendextend)
library(corrplot)
library(flexclust)

# Za³adowanie biblioteki dla analityki tematyk
library(topicmodels)

# Za³adowanie biblioteki dla analizy s³ów kluczowych
library(wordcloud)

# Zmiana katalogu roboczego
workDir <- "D:\\Software\\UEK\\PJN_Projekt"
setwd(workDir)

# Definicja œcie¿ek dostêpu do katalogów funkcjonalnych
inputDir <- ".\\Dane"
outputDir <- ".\\Wyniki"

# Przetwarzanie wstêpne (podpunkt A oraz B) ---------------------------------------------------


# Utworzenie korpusu dokumentów dla przetwarzania wstêpnego
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

# Wstêpne przetwarzanie
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

# Usuniêcie znaków em dash i 3/4
removeCharPP <- content_transformer(
    function(x,character) gsub(character, "",x)
)
corpusPP <- tm_map(corpusPP, removeCharPP, intToUtf8(8722))
corpusPP <- tm_map(corpusPP, removeCharPP, intToUtf8(190))

# Usuniêcie rozszerzeñ z nazw dokumentów w korpusie
cutExtensionsPP <- function(document){
    meta(document,"id") <- gsub(
        pattern = "\\.txt$",
        replacement = "",
        meta(document,"id")
    )
    return(document)
}
corpusPP <- tm_map(corpusPP, cutExtensionsPP)

#usuniêcie podzia³u na akapity w dokumentach tekstowych
pasteParagraphsPP <- content_transformer(
    function(text, char) paste(text, collapse = char)
)
corpusPP <- tm_map(corpusPP, pasteParagraphsPP, " ")

#lematyzacja
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

#eksport przetworzonego korpusu
preprocessedDir <- paste(
    outputDir,
    "Dokumenty_Przetworzone",
    sep = "\\"
)
dir.create(preprocessedDir, showWarnings = FALSE)
writeCorpus(corpusPP, path = preprocessedDir)


# Macierze czêstoœci (podpunkt C) ------------------------------------------------------


# Utworzenie korpusu dokumentów
corpusMatrixDir <- paste(
    outputDir,
    "Dokumenty_Przetworzone",
    sep = "\\"
)

corpusMatrix <- VCorpus(
    DirSource(
        corpusMatrixDir,
        "CP1250",
        "*.txt"
    ),
    readerControl = list(
        language = "pl_PL"
    )
)

# Usuniêcie rozszerzeñ z nazw dokumentów w korpusie
cutExtensionsMatrix <- function(document){
    meta(document,"id") <- gsub(
        pattern = "\\.txt$",
        replacement = "",
        meta(document,"id")
    )
    return(document)
}
corpusMatrix <- tm_map(corpusMatrix, cutExtensionsMatrix)

# Tworzenie macierzy czêstoœci
tdmTfAll <- TermDocumentMatrix(corpusMatrix)
tdmTfIdfAll <- TermDocumentMatrix(
    corpusMatrix,
    control = list(
        weighting = weightTfIdf
    )
)
tdmTfBounds <- TermDocumentMatrix(
    corpusMatrix,
    control = list(
        bounds = list(
            global = c(2,16)
        )
    )
)
tdmTfIdfBounds <- TermDocumentMatrix(
    corpusMatrix,
    control = list(
        weighting = weightTfIdf,
        bounds = list(
            global = c(2,16)
        )
    )
)

dtmTfAll <- DocumentTermMatrix(corpusMatrix)
dtmTfIdfAll <- DocumentTermMatrix(
    corpusMatrix,
    control = list(
        weighting = weightTfIdf
    )
)
dtmTfBounds <- DocumentTermMatrix(
    corpusMatrix,
    control = list(
        bounds = list(
            global = c(2,16)
        )
    )
)
dtmTfIdfBounds <- DocumentTermMatrix(
    corpusMatrix,
    control = list(
        weighting = weightTfIdf,
        bounds = list(
            global = c(2,16)
        )
    )
)

# Konwersje macierzy rzadkich do macierzy klasycznych
tdmTfAllMatrix <- as.matrix(tdmTfAll)
tdmTfIdfAllMatrix <- as.matrix(tdmTfIdfAll)
tdmTfBoundsMatrix <- as.matrix(tdmTfBounds)
tdmTfIdfBoundsMatrix <- as.matrix(tdmTfIdfBounds)
dtmTfAllMatrix <- as.matrix(dtmTfAll)
dtmTfIdfAllMatrix <- as.matrix(dtmTfIdfAll)
dtmTfBoundsMatrix <- as.matrix(dtmTfBounds)
dtmTfIdfBoundsMatrix <- as.matrix(dtmTfIdfBounds)

# Eksport macierzy czêstoœci do pliku
matrixFile <- paste(
  outputDir,
  "tdmTfAllMatrix.csv",
  sep = "\\"
)
write.table(
  tdmTfAllMatrix, 
  file = matrixFile, 
  sep = ";", 
  dec = ",", 
  col.names = NA
)


# Redukcje wymiarów (podpunkt D) -------------------------------------------------------


##-- Redukcja przy u¿yciu analizy g³ównych sk³adowych --##

# Analiza g³ównych sk³adowych dla wczeœniej utworzonej macierzy czêstoœci
pca <- prcomp(dtmTfIdfBounds)

#przygotowanie danych do wykresu
legendPCA <- paste(paste("d",1:20,sep = ""), rownames(dtmTfIdfBounds), sep = " => ")
xPCA <- pca$x[,1]
yPCA <- pca$x[,2]

#eksport wykresu do pliku
plotFilePCA <- paste(
    outputDir,
    "pca.png",
    sep = "\\"
)
png(filename = plotFilePCA, width = 1920, height = 1080)

#wykres w przestrzeni dwuwymiarowej
plot(
    xPCA,
    yPCA, 
    main = "Analiza g³ównych sk³adowych",
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

#zamkniêcie pliku
dev.off()

##-- Redukcja przy u¿yciu dekompozycji wed³ug wartoœci osobliwych --##

# Analiza ukrytych wymiarów semantycznych 
#(dekompozycja wg wartoœci osobliwych)

lsa <- lsa(tdmTfAllMatrix)

#przygotowanie danych do wykresu
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

#eksport wykresu do pliku
plotFileLSA <- paste(
    outputDir,
    "lsa.png",
    sep = "\\"
)
png(filename = plotFileLSA, height = 1920, width = 1080)

#wykres w przestrzeni dwuwymiarowej
plot(
    x1LSA,
    y1LSA, 
    main = "Analiza ukrytych wymiarów semantycznych",
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

#zamkniêcie pliku
dev.off()


# Analiza skupieñ dokumentów (podpunkt E) ----------------------------------------------

##-- Metoda hierarchiczna --##
##-- Parametry metody --##
# 1. Macierz czêstoœci
# a. Waga (weighting)
# b. Z zakresem uwzglêdnionych zmiennych (bounds)
# 2. Miara odleg³oœci (euclidean, jaccard, cosine)
# 3. Sposób wyznaczania odleg³oœci pomiêdzy skupieniami 
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
nTermsLDA <- ncol(dtmTfAll)
nTopicsLDA <- 4
lda <- LDA(
    dtmTfAll,
    k = nTopicsLDA,
    method = "Gibbs",
    control = list(
        burnin = 2000,
        thin = 100,
        iter = 3000
    )
)
resultsLDA <- posterior(lda)

# Prezentacja tematów
par(mai = c(1,2,1,1))
coloursLDA <- c("violet","orange","turquoise","darkseagreen")
topic1 <- head(sort(resultsLDA$terms[1,],decreasing = TRUE),20)
barplot(
    rev(topic1),
    horiz = TRUE,
    las = 1, 
    main = "Temat 1",
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA[1]
)
topic2 <- head(sort(resultsLDA$terms[2,],decreasing = TRUE),20)
barplot(
    rev(topic2),
    horiz = TRUE,
    las = 1, 
    main = "Temat 2",
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA[2]
)
topic3 <- head(sort(resultsLDA$terms[3,],decreasing = TRUE),20)
barplot(
    rev(topic3),
    horiz = TRUE,
    las = 1, 
    main = "Temat 3",
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA[3]
)
topic4 <- head(sort(resultsLDA$terms[4,],decreasing = TRUE),20)
barplot(
    rev(topic4),
    horiz = TRUE,
    las = 1, 
    main = "Temat 4",
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA[4]
)

# Prezentacja dokumentów
document1 <- resultsLDA$topics[1,]
barplot(
    document1,
    horiz = TRUE,
    las = 1, 
    main = rownames(resultsLDA$topics)[1],
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA
    #col = coloursLDA[which.max(resultsLDA$topics[1,])]
)
document4 <- resultsLDA$topics[4,]
barplot(
    document4,
    horiz = TRUE,
    las = 1, 
    main = rownames(resultsLDA$topics)[4],
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA
    #col = coloursLDA[which.max(resultsLDA$topics[4,])]
)
document11 <- resultsLDA$topics[11,]
barplot(
    document11,
    horiz = TRUE,
    las = 1, 
    main = rownames(resultsLDA$topics)[11],
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA
    #col = coloursLDA[which.max(resultsLDA$topics[11,])]
)
document17 <- resultsLDA$topics[17,]
barplot(
    document17,
    horiz = TRUE,
    las = 1, 
    main = rownames(resultsLDA$topics)[17],
    xlab = "Prawdopodobieñstwo",
    col = coloursLDA
    #col = coloursLDA[which.max(resultsLDA$topics[17,])]
)

# Udzia³ tematów w s³owach
options(scipen = 5)
words1 <- c("pandemia")
round(resultsLDA$terms[,words1],4)

words2 <- c("ekonomia")
round(resultsLDA$terms[,words2],4)

words3 <- c("zdrowie")
round(resultsLDA$terms[,words3],4)


# Analiza s³ów kluczowych (podpunkt G) ----------------------------------------------------------

#-- Dla pierwszego dokumentu --#
#-- Waga tf jako miara wa¿noœci s³ów --#
keywordsTf1 <- head(sort(dtmTfBoundsMatrix[1,],decreasing = TRUE))
keywordsTf1

# -- Waga tfidf jako miara wa¿noœci s³ów --#
keywordsTfIdf1 <- head(sort(dtmTfIdfBoundsMatrix[1,],decreasing = TRUE))
keywordsTfIdf1

#-- Prawdopodobieñstwo w modelu LDA jako miara wa¿noœci s³ów --##
termsImportance1 <- c(resultsLDA$topics[1,]%*%resultsLDA$terms)
names(termsImportance1) <- colnames(resultsLDA$terms)
keywordsLda1 <- head(sort(termsImportance1,decreasing = TRUE))
keywordsLda1

#-- Chmura tagów --##
par(mai = c(0,0,0,0))
wordcloud(corpusPP[2], max.words = 200, colors = brewer.pal(8,"PuOr"))