library(Mfuzz)
library(Biobase)

dataset_1 <- read.csv("E:/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua proteomics/Clusters/Clusters final/dataset_1_ok.csv",sep=";", row.names=1)
dataset_2 <- read.csv("E:/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua proteomics/Clusters/Clusters final/dataset_2.csv",sep=";", row.names=1)






eset <- ExpressionSet(assayData = as.matrix(dataset_1),
                      phenoData = AnnotatedDataFrame(dataset_2))

str(eset)


eset <- standardise(eset)



cl12 <- mfuzz(eset,c=7,m=2.00)
mfuzz.plot(eset,cl=cl12,mfrow=c(4,4))
O <- overlap(cl) 
Ptmp <- overlap.plot(cl,over=O,thres=0.05)

O3 <- overlap(cl12)
overlap.plot(cl12,over=O3,thres=0.05,P=Ptmp)
cl[[4]][1:20,1]
cl[[4]][1:20,4]
cl3 <- mfuzz(eset =cl,c=10,m=3.85)

acore(eset, cl=cl12, min.acore = 0)

m1 <- mestimate(eset)
m1

help ("mfuzz")


write.csv(oLINK_1,"C:/Users/USUARIO/Desktop/average_sig.csv")


Mfuzzgui()
