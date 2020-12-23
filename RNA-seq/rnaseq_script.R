library(sva)
library(pamr)
library(limma)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#set directory with input data
homedir = "immuno"
plot_dir = "plots"

if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
} else {
  print("Dir already exists!")
}

#read gene counts
counts = read.table(paste0(homedir,'/featurecounts.txt'), header=T, row.names = 1)

#read phenotype file
pheno = read.table(paste0(homedir,'/pheno.csv'), header=T, sep=',', row.names = 1)


#normalize counts to log2(CPM)
design.block = model.matrix(~celltype, data=pheno)
voomExp = voom(counts, design=design.block, plot=T)

#create model for ComBat taking into account cell type and run normalization
combatmod = model.matrix(~1+celltype, data=pheno)

edata_combat = ComBat(dat=voomExp$E, batch=pheno$batch, mod=combatmod, par.prior = T, prior.plots = T)

edata_combat = as.data.frame(edata_combat)

#read annotation file with gene IDs and names
selected_genes = read.table(paste0(homedir,'/genes.txt'),header=T)

#subset gene IDs and names
ngene <- selected_genes[,9:10]

selected_counts = subset(edata_combat, row.names(edata_combat) %in% ngene$gene_id)

selected_counts$Gene <- rownames(selected_counts)


t <- rownames(pheno[pheno$celltype=="T-cell",])
b <- rownames(pheno[pheno$celltype=="B-cell",])
nk<- rownames(pheno[pheno$celltype=="NK",])

genes <- as.vector(selected_counts$Gene)

out = list()
  
#prepare plots
for (g in genes){
  #subset expression for selected gene into new gen1 data frame
  gen1 <- data.frame(selected_counts[selected_counts$Gene==g,],row.names = NULL )
  gname <- as.character(ngene[as.character(ngene$gene_id)==g,]$gene_name)
  
  gg <- data.frame(t(gen1))
  gg$samples <- rownames(gg)
  rownames(gg) <- NULL
  gg <- gg[1:nrow(gg)-1,]
  
  #add group column to gen1 data frame
  gg$group[gg$samples %in% b] <- "B-cells"
  gg$group[gg$samples %in% t] <- "T-cells"
  gg$group[gg$samples %in% nk] <- "NK-cells"
  
  gg$t.gen1. = as.numeric(as.character(gg$t.gen1.))
  gg$group = as.factor(gg$group)
  
  #check if all rows have group label
  gg$group <- factor(gg$group , levels=c("T-cells", "B-cells", "NK-cells"))
  gg = gg[gg$group %in% c("B-cells","NK-cells","T-cells") ,]
  
  #fit an analysis of variance model
  comps = list(c("T-cells", "B-cells"),c("T-cells","NK-cells"),c("B-cells","NK-cells"))
  a = aov(gg$t.gen1. ~ gg$group)
  aov_p = summary(a)[[1]][["Pr(>F)"]]
  
  range <- max(gg$t.gen1.)-min(gg$t.gen1.)
  
  #plot boxplots
  p <- ggplot(gg, aes(x=group, y=t.gen1., fill=group))  + 
    geom_boxplot(outlier.color="red", outlier.shape=4, outlier.size=4)+ theme_classic() +
    stat_compare_means(step.increase=0.2, label = "p.signif", method = "t.test",comparisons=comps) + 
    scale_x_discrete(labels=c("T-cells" = "T cells", "B-cells" = "B cells", "NK-cells" = "NK cells")) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=4) + ylab(paste0(gname," mRNA level [log2]\n"))
    
  png(file=paste(plot_dir,"/",gname,".png",sep=""))
  print(p)
  dev.off()
  
}


#plot PRDX1 expression for selected 4 samples 

ngene[ngene$gene_name=="PRDX1",]
sample1 <- c("SRR1550988","SRR1550989","SRR1550990","SRR1550991")
sample2 <- c("SRR1551049","SRR1551050","SRR1551051","SRR1551052")
sample3 <- c("SRR1551070","SRR1551071","SRR1551072","SRR1551073")
sample4 <- c("SRR1551056","SRR1551057","SRR1551058","SRR1551059")

prdx1 <- selected_counts[selected_counts$Gene==ngene[ngene$gene_name=="PRDX1",'gene_id'],]


sample1_counts <- prdx1[,colnames(prdx1) %in% sample1]
sample2_counts <- prdx1[,colnames(prdx1) %in% sample2]
sample3_counts <- prdx1[,colnames(prdx1) %in% sample3]
sample4_counts <- prdx1[,colnames(prdx1) %in% sample4]

#for T-cells get average from CD4 and CD8 cells
sample1_counts$T <- (sample1_counts$SRR1550989+sample1_counts$SRR1550990)/2
sample2_counts$T <- (sample2_counts$SRR1551050+sample2_counts$SRR1551051)/2
sample3_counts$T <- (sample3_counts$SRR1551071+sample3_counts$SRR1551072)/2
sample4_counts$T <- (sample4_counts$SRR1551057+sample4_counts$SRR1551058)/2


sample1_counts <- sample1_counts[,c(1,4:5)]
colnames(sample1_counts) <- c("B cells", "NK cells", "T cells")


sample2_counts <- sample2_counts[,c(1,4:5)]
colnames(sample2_counts) <- c("B cells", "NK cells", "T cells")


sample3_counts <- sample3_counts[,c(1,4:5)]
colnames(sample3_counts) <- c("B cells", "NK cells", "T cells")


sample4_counts <- sample4_counts[,c(1,4:5)]
colnames(sample4_counts) <- c("B cells", "NK cells", "T cells")


prdx1_final <- rbind(sample1_counts,sample2_counts,sample3_counts,sample4_counts)

patients <- c("A","B","C","D")
prdx1_final$patients <- patients

#create final data frame with expression
total_exp <- c(prdx1_final$`T cells`,prdx1_final$`B cells`,prdx1_final$`NK cells`)
cells <- c(rep("T cells",4),rep("B cells",4), rep("NK cells",4))
pats <- c(rep(c("A","B","C","D"),3))

final_df <- data.frame(prdx1=total_exp,cells=cells,pats=pats)

final_df$cells <- factor(final_df$cells , levels=c("T cells", "B cells", "NK cells"))

q <- ggplot(final_df, aes(x=cells, y=prdx1, color=pats, group=pats))+geom_point(size=10) + geom_line(size=1,linetype="dotted") +
  theme_classic() + ylab("PRDX1 mRNA level [log2]\n")

png(file=paste(plot_dir,"/PRDX1_dots",".png",sep=""), height = 700, width = 700)
print(q)
dev.off()
