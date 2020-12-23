library(edgeR)
library(ggpubr)


#set directory with input data
homedir = "immuno/mapped/"
plot_dir = "plots_nk"

if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
} else {
  print("Dir already exists!")
}


#get list of gene count files
ff = list.files(path=homedir, pattern = ".*.txt")

#combine raw counts for all samples
t = read.table(paste0(homedir,ff[1]), header=F, row.names = 1, sep='\t', fill=TRUE)
t <- t[1:(nrow(t)-5),]
gdup <- t$V2[duplicated(t$V2)]
t <- t[!(t$V2 %in% gdup),]

countsRaw <- data.frame(row.names = t$V2)

for (f in ff) {
  tmp = read.table(paste0("mapped/",f), header=F, row.names = 1, sep='\t', fill=TRUE)
  tmp <- tmp[1:(nrow(tmp)-5),]
  tmp <- tmp[tmp$V2 %in% rownames(countsRaw),]
  countsRaw <- cbind(countsRaw, tmp$V3)
  print(f)
}

#set sample names
samples <- c("IL15_W_3B","IL2_W_1","IL15_W_1","IL15_W_2","IL2_A_2","IL15_A_4","IL2_A_4","IL2_A_1","IL2_W_4","IL15_A_5","IL15_W_5","IL2_A_3A","IL2_W_2","IL2_W_5","IL15_A_3A",
             "IL2_W_3B","IL2_A_5","IL15_W_4","IL15_A_2","IL15_A_1")

colnames(countsRaw) <- samples


#get vector with sample groups
condition <- c()
for (i in colnames(cpms)){
  out = strsplit(i,"_")[[1]]
  res <- paste0(out[1],"_",out[2])
  condition <- c(condition, res)
}



#use TMM normalization 
y = DGEList(counts=countsRaw)
y.norm = calcNormFactors(y, method="TMM")

#normalize expression to CPM
cpms = cpm(y.norm, group=condition, log=FALSE, normalized.lib.sizes=TRUE)


#select genes for analysis
outputGenes <- c("PRDX1", "TXN", "TXNRD1","TXNIP")

#select groups to compare
ind = (condition=="IL15_A" | condition=="IL15_W")

#get subset of data for analysis
cpmsG <- cpms[,ind]
groups <- condition[ind]
cpmsG <- cpmsG[outputGenes,]


#prepare plots
for (gen in rownames(cpmsG)){
  #subset expression for selected gene to new data frame dat
  vals <- as.data.frame(t(cpmsG))[,gen]
  dat <- data.frame(vals = vals, group = groups)
  comps = list(c("IL15_A", "IL15_W"))
  p <- ggplot(dat, aes(x=group, y=log2(vals+1), fill=group))  + 
    geom_boxplot(outlier.color="red", outlier.shape=4, outlier.size=4)+ theme_classic() +
    theme_classic() +  ylab(paste0(gen," polysome level\n[log2(CPM)]\n")) +
    stat_compare_means(label = "p.signif", method = "wilcox.test", paired=TRUE, comparisons=comps) 
  png(file=paste0(plot_dir,"/",gen,".png"))
  print(p)
  dev.off()
  
}



