% Biochem 158 Final Project (TITLE)
% Kris Sankaran
% June 7, 2013

Data Analysis
--------------

Here, we reproduce analysis of the paper Witten etc.






```r
miRNA.counts <- read.csv("nSequencesObtained.csv", header = TRUE)
head(miRNA.counts)  # These are the number of reads in each individual
```

```
  SampleID Diagnosis Normal  Cancer
1     G547       ADC 335672 1329659
2     G659       ADC 911955  399605
3     G691       ADC 625927  185932
4     G696       ADC 254493   29848
5     G761       ADC 638140   71636
6     G220       ADC  65811  391371
```

```r
library(ggplot2)

library(reshape2)
mCounts <- melt(miRNA.counts[, c("Normal", "Cancer")], variable.name = "type", 
    value.name = "nReads")
```

```
Using as id variables
```

```r
sort(mCounts[, "nReads"], decreasing = TRUE)
```

```
 [1] 25007613 17340713  2624426  1909053  1796567  1737419  1707320
 [8]  1650764  1404909  1340962  1329659  1218235  1192819  1174576
[15]  1137870  1072585  1012212   971428   949415   937452   911955
[22]   820035   756800   750034   745293   742291   736372   679950
[29]   661193   638140   632680   625927   599073   578876   508609
[36]   491186   463677   399605   391371   389491   367475   339206
[43]   335672   319933   308769   304862   254493   219426   185932
[50]   180675   160869   133790   125175    91461    90111    74780
[57]    71636    65811    59984    29848
```

```r
ggplot(mCounts[-c(30, 45, 60), ]) + geom_bar(aes(x = nReads, fill = type), binwidth = 1e+05) + 
    ggtitle("Histogram of Counts per Sample")  # Excluding largest counts, for viewability
```

![plot of chunk overallSequencesCounts](figure/overallSequencesCounts.png) 



```r
miRNA.expression <- read.csv("miRNA_expression.csv", header = TRUE)
mExpress <- melt(miRNA.expression, variable.name = "type", value.name = "nReads")
```

```
Using miRNA as id variables
```

```r
mExpress$nReads <- (mExpress$nReads)^(1/3)  # cube root transformation
normal.ix <- which(substr(mExpress[, "type"], 1, 1) == "N")
mExpressSummary <- mExpress
mExpressSummary[normal.ix, "type"] <- "N"
mExpressSummary[-normal.ix, "type"] <- "T"

ggplot(mExpressSummary) + geom_histogram(aes(x = nReads, fill = type)) + ggtitle("Number of Reads in Tumor and Normal Samples")
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

![plot of chunk countsPerSampleType](figure/countsPerSampleType.png) 



```r
X <- read.csv("miExpress2.csv", header = T)
rownames(X) <- X[, 1]
X <- X[, -1]

X <- X^(1/3)  # Cube-root transform
miRNA.clust <- hclust(as.dist(1 - cor(X)))
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")  # Plotting Package
A2Rplot(miRNA.clust, k = 3, col.down = c("#FF6B6B", "#4ECDC4", "#556270"), main = "Hierchical Clustering of Patients using miRNA Counts")
```

![plot of chunk hClust](figure/hClust.png) 



```r
X.tilde <- scale(t(X), center = TRUE, scale = FALSE)
svdX <- svd(X.tilde)
pca.X <- svdX$u %*% diag(svdX$d)
pca.approx <- pca.X[, 1:2]

## Plot PCA
colnames(pca.approx) <- c("PC1", "PC2")
pca.labeled <- data.frame(pca.approx, label = colnames(X))
sample.type <- substr(pca.labeled$label, start = 6, stop = 6)
sample.type[59:60] <- c("N", "T")
ggplot(pca.labeled) + geom_text(aes(x = PC1, y = PC2, label = label, col = sample.type)) + 
    ggtitle("Principal Components from miRNAs Discriminate between Groups")
```

![plot of chunk pca](figure/pca.png) 

