

#########################
# Starting with the c7atoms and the c7 ssGSEA scores, getting the top correlates out.
########################

genesets <- read.table("msigdb_c7_immune_gene_sets.txt")
humansets <- read.table("human_c7_genesets.txt")
genesets <- genesets[genesets$V1 %in% humansets$V1,]

load("data/Signatures/ssGSEA_Scores_Wolf68_Bindea_Yasin_C7.rda") # faster
load("data/feature_starting_points/Signature_Components.rda")

c7data <- res0[res0$Source=="MSigDB_C7", -c(1,2)]

rownames(c7data) <- res0[res0$Source=="MSigDB_C7", 2]

c7dt <- t(c7data)

c7at <- t(c7atoms)

m <- cor(c7at, c7dt, method="spearman")

top3cor <- lapply(1:nrow(m), function(i) {m[i, which( abs(m[i,]) >= sort(abs(m[i,]),decreasing=T)[3])]} )

top3names <- lapply(1:nrow(m), function(i) {colnames(m)[which( abs(m[i,]) >= sort(abs(m[i,]),decreasing=T)[3])]} )


genesets[genesets$V1 == "GSE34515_CD16_POS_MONOCYTE_VS_DC_DN",2]
write.table(genesets[genesets$V1 == "GSE34515_CD16_POS_MONOCYTE_VS_DC_DN",2], file="gse34515_genes.txt", quote=F, row.names=F, col.names=F)
