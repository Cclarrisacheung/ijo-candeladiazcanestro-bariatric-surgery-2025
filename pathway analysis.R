if (!require("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")
BiocManager::install(version = "3.16")


df = read.csv("E:/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua proteomics/Pathway analysis/1M.wt.csv",header = TRUE, sep = ";")

original_gene_list <- df$log2FoldChange


names(original_gene_list) <- df$Symbol


gene_list<-na.omit(original_gene_list)


gene_list = sort(gene_list, decreasing = TRUE)

gene_list


library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(pathview)
library(ReactomePA)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

gse

require(DOSE)
dotplot(gse, showCategory=30, split=".sign") + facet_grid(.~.sign)
emapplot(gse, showCategory = 30)

write.csv(ids ,"C:/Users/USUARIO/Desktop/ids.csv")


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(df$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db")
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$Symbol %in% dedup_ids$SYMBOL,]
df2
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y
names(kegg_gene_list)
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
emapplot(kk2)
library(enrichplot)
barplot(kk2, showCategory=20)
cnetplot(kk2, categorySize="geneNum", foldChange=gene_list, showCategory = 20)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list,15)

library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism)
