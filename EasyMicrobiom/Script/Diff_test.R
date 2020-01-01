
# DESep2(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")

DESep2 = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv"){
  library(phyloseq)
  #导入otu表格
  otu = read.delim(otu,row.names = 1)
  head(otu)
  otu = as.matrix(otu)
  str(otu)
  #导入注释文件
  tax = read.delim(tax,row.names = 1)
  head(tax)
  tax = as.matrix(tax)
  # taxa_names(tax)
  
  #导入分组文件
  map = read.delim(map,row.names = 1)
  head(map)
  
  # #导入进化树
  # tree = read.tree("./otus.tree")
  # tree
  
  ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                 sample_data(map) ,
                 tax_table(tax)
                 # phy_tree(tree)
                 
  )
  ps
  
  
  sub_design <- as.data.frame(sample_data(ps))
  Desep_group <- as.character(levels(sub_design$Group))
  Desep_group
  
  aaa = combn(Desep_group,2)
  ##设置文件路径
  path = "./DESep2/"
  dir.create(path)
  # aaa[,1]
  # i = 1
  for (i in 1:dim(aaa)[2]) {
    
    Desep_group = aaa[,i]
    
    library("DESeq2")
    pkgs <- c("phyloseq", "structSSI", "dplyr", "reshape2",
              "ggplot2", "DESeq2")
    sapply(pkgs, require, character = TRUE)
    vegan_otu <-  function(physeq){
      OTU <-  otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU <-  t(OTU)
      }
      return(as(OTU,"matrix"))
    }
    otu_table = as.data.frame(vegan_otu(ps))
    
    count = as.matrix(otu_table)
    count <- t(count)
    #数据整理形式一######
    
    sub_design <- as.data.frame(sample_data(ps))
    dim(sub_design)
    sub_design$SampleType = as.character(sub_design$Group)
    sub_design$SampleType <- as.factor(sub_design$Group)
    
    
    
    
    dds <- DESeqDataSetFromMatrix(countData = count,
                                  colData = sub_design,
                                  design = ~ SampleType)
    
    dds2 <- DESeq(dds)  ##第二步,标准化
    resultsNames(dds2)
    # 将结果用results()函数来获取，赋值给res变量
    res <-  results(dds2, contrast=c("SampleType",Desep_group ),alpha=0.05)
    # summary一下，看一下结果的概要信息
    summary(res)
    # 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
    WT <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
    
    dim(WT)
    res$level = as.factor(ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "enriched",ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "depleted","nosig")))
    # dim(res)
    # head(res)
    x <- res
    ###########添加物种丰度
    # dim(count)
    # str(count)
    count = as.matrix(count)
    norm = t(t(count)/colSums(count)) #* 100 # normalization to total 100
    dim(norm)
    norm1 = norm %>% 
      t() %>% as.data.frame()
    # head(norm1)
    #数据分组计算平均值
    library("tidyverse")
    dim(norm1)
    
    iris.split <- split(norm1,as.factor(sub_design$SampleType))
    iris.apply <- lapply(iris.split,function(x)colMeans(x))
    # 组合结果
    iris.combine <- do.call(rbind,iris.apply)
    norm2= t(iris.combine)
    
    #head(norm)
    str(norm2)
    norm2 = as.data.frame(norm2)
    # dim(x)
    # dim(norm2)
    x = cbind(x,norm2)
    # head(x)
    #在加入这个文件taxonomy时，去除后面两列不相干的列
    # 读取taxonomy，并添加各列名称
    vegan_tax <-  function(physeq){
      tax <-  tax_table(physeq)
      
      return(as(tax,"matrix"))
    }
    taxonomy = as.data.frame(vegan_tax(ps))
    head(taxonomy)
    # taxonomy <- as.data.frame(tax_table(ps1))
    
    #发现这个注释文件并不适用于直接作图。
    #采用excel将其分列处理，并且删去最后一列，才可以运行
    if (length(colnames(taxonomy)) == 6) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    }else if (length(colnames(taxonomy)) == 7) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
    }else if (length(colnames(taxonomy)) == 8) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","rep")
    }
    # colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    
    # Taxonomy排序，并筛选OTU表中存在的
    library(dplyr)
    taxonomy$id=rownames(taxonomy)
    # head(taxonomy)
    tax = taxonomy[row.names(x),]
    
    if (length(colnames(taxonomy)) == 7) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      # x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 8) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 9) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }
    
    head(x)
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    
    filename = paste(path,"/a3_Desep_for_",group,"_allOTU_.csv",sep = "")
    write.csv(x,filename,quote = FALSE)
    
    filename = paste(path,"/a3_Desep_for_",group,"differ_OTU.csv",sep = "")
    write.csv(WT,filename,quote = FALSE)
    
    
  }
  
}

# Edger(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")

Edger = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv"){
  library(phyloseq)
  #导入otu表格
  otu = read.delim(otu,row.names = 1)
  head(otu)
  otu = as.matrix(otu)
  str(otu)
  #导入注释文件
  tax = read.delim(tax,row.names = 1)
  head(tax)
  tax = as.matrix(tax)
  # taxa_names(tax)
  
  #导入分组文件
  map = read.delim(map,row.names = 1)
  head(map)
  
  # #导入进化树
  # tree = read.tree("./otus.tree")
  # tree
  
  ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                 sample_data(map) ,
                 tax_table(tax)
                 # phy_tree(tree)
                 
  )
  ps
  sub_design <- as.data.frame(sample_data(ps))
  Desep_group <- as.character(levels(sub_design$Group))
  Desep_group
  aaa = combn(Desep_group,2)
  sub_design <- as.data.frame(sample_data(ps))
  path = "./Edgr"
  
  dir.create(path)
  # aaa[,1]
  for (i in 1:dim(aaa)[2]) {
    
    Desep_group = aaa[,i]
    
    library("DESeq2")
    pkgs <- c("phyloseq", "structSSI", "dplyr", "reshape2",
              "ggplot2", "DESeq2")
    sapply(pkgs, require, character = TRUE)
    vegan_otu <-  function(physeq){
      OTU <-  otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU <-  t(OTU)
      }
      return(as(OTU,"matrix"))
    }
    otu_table = as.data.frame(vegan_otu(ps))
    
    count = as.matrix(otu_table)
    count <- t(count)
    #数据整理形式一######
    
    sub_design <- as.data.frame(sample_data(ps))
    dim(sub_design)
    sub_design$SampleType = as.character(sub_design$Group)
    sub_design$SampleType <- as.factor(sub_design$Group)
    library(edgeR)
    
    # create DGE list
    d = DGEList(counts=count, group=sub_design$SampleType)
    d$samples
    d = calcNormFactors(d)#默认为TMM标准化
    # 生成实验设计矩阵
    design.mat = model.matrix(~ 0 + d$samples$group)
    colnames(design.mat)=levels(sub_design$SampleType)
    d2 = estimateGLMCommonDisp(d, design.mat)
    d2 = estimateGLMTagwiseDisp(d2, design.mat)
    fit = glmFit(d2, design.mat)
    # head(design)
    # 设置比较组写在前面的分组为enrich表明第一个分组含量高
    
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    BvsA <- makeContrasts(contrasts =  group,levels=design.mat)#注意是以GF1为对照做的比较
    # 组间比较,统计Fold change, Pvalue
    lrt = glmLRT(fit,contrast=BvsA)
    
    # FDR检验，控制假阳性率小于5%
    de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05,lfc=0)#lfc=0这个是默认值
    summary(de_lrt)
    # 导出计算结果
    x=lrt$table
    x$sig=de_lrt
    # head(x)
    x <- cbind(x, padj = p.adjust(x$PValue, method = "fdr"))
    enriched = row.names(subset(x,sig==1))
    depleted = row.names(subset(x,sig==-1))
    x$level = as.factor(ifelse(x$sig==1, "enriched",ifelse(x$sig==-1, "depleted","nosig")))
    head(x)
    WT <-subset(x,padj < 0.05 )
    dim(WT)
    
    
    
    ###########添加物种丰度
    # dim(count)
    # str(count)
    count = as.matrix(count)
    norm = t(t(count)/colSums(count)) #* 100 # normalization to total 100
    dim(norm)
    norm1 = norm %>% 
      t() %>% as.data.frame()
    # head(norm1)
    #数据分组计算平均值
    library("tidyverse")
    dim(norm1)
    
    iris.split <- split(norm1,as.factor(sub_design$SampleType))
    iris.apply <- lapply(iris.split,function(x)colMeans(x))
    # 组合结果
    iris.combine <- do.call(rbind,iris.apply)
    norm2= t(iris.combine)
    
    #head(norm)
    str(norm2)
    norm2 = as.data.frame(norm2)
    # dim(x)
    # dim(norm2)
    x = cbind(x,norm2)
    # head(x)
    #在加入这个文件taxonomy时，去除后面两列不相干的列
    # 读取taxonomy，并添加各列名称
    vegan_tax <-  function(physeq){
      tax <-  tax_table(physeq)
      
      return(as(tax,"matrix"))
    }
    taxonomy = as.data.frame(vegan_tax(ps))
    head(taxonomy)
    # taxonomy <- as.data.frame(tax_table(ps1))
    
    #发现这个注释文件并不适用于直接作图。
    #采用excel将其分列处理，并且删去最后一列，才可以运行
    if (length(colnames(taxonomy)) == 6) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    }else if (length(colnames(taxonomy)) == 7) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
    }else if (length(colnames(taxonomy)) == 8) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","rep")
    }
    # colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    
    # Taxonomy排序，并筛选OTU表中存在的
    library(dplyr)
    taxonomy$id=rownames(taxonomy)
    # head(taxonomy)
    tax = taxonomy[row.names(x),]
    
    if (length(colnames(taxonomy)) == 7) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      # x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 8) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 9) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }
    
    head(x)
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    
    filename = paste(path,"/a3_Desep_for_",group,"_allOTU_.csv",sep = "")
    write.csv(x,filename,quote = FALSE)
    
    filename = paste(path,"/a3_Desep_for_",group,"differ_OTU.csv",sep = "")
    write.csv(WT,filename,quote = FALSE)
    
    
  }
}


# Ttest(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
Ttest = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv"){
  #导入otu表格
  library(phyloseq)
  otu = read.delim(otu,row.names = 1)
  head(otu)
  otu = as.matrix(otu)
  str(otu)
  #导入注释文件
  tax = read.delim(tax,row.names = 1)
  head(tax)
  tax = as.matrix(tax)
  # taxa_names(tax)
  
  #导入分组文件
  map = read.delim(map,row.names = 1)
  head(map)
  
  # #导入进化树
  # tree = read.tree("./otus.tree")
  # tree
  
  ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                 sample_data(map) ,
                 tax_table(tax)
                 # phy_tree(tree)
                 
  )
  ps
  
  ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  sub_design <- as.data.frame(sample_data(ps))
  
  
  Desep_group <- as.character(levels(sub_design$Group))
  Desep_group
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  otu_table = as.data.frame(vegan_otu(ps_rela))
  
  count = as.matrix(otu_table)
  count <- t(count)
  #数据整理形式一######
  
  sub_design <- as.data.frame(sample_data(ps_rela))
  dim(sub_design)
  sub_design$SampleType = as.character(sub_design$Group)
  sub_design$SampleType <- as.factor(sub_design$Group)
  
  aaa = combn(Desep_group,2)
  ##设置文件路径
  path = "./T_test"
  
  dir.create(path)
  # aaa[,1]
  # i =1
  for (i in 1:dim(aaa)[2]) {
    Desep_group = aaa[,i]
    # 转换原始数据为百分比，
    norm = t(t(count)/colSums(count,na=TRUE)) * 100 # normalization to total 100
    head(norm)
    norm=as.data.frame(norm)
    a=norm
    head(a)
    #预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
    Pvalue<-c(rep(0,nrow(a)))
    fdr<-c(rep(0,nrow(a)))
    log2_FC<-c(rep(0,nrow(a)))
    sub_design$ID = row.names(sub_design)
    library(tidyverse)
    df_filter<- filter(sub_design ,Group %in% Desep_group)
    
    head(df_filter)
    
    a = as.data.frame(a)
    a = a[as.character(df_filter$ID)]
    head(a)
    rep = length(as.character(df_filter$ID))/2
    
    a = as.matrix(a)
    ###########开始运行脚本
    for(i in 1:nrow(a)){
      if(sd(a[i,(1:rep)])==0&&sd(a[i,(rep+1):(rep*2)])==0){
        Pvalue[i] <-"NA"
        log2_FC[i]<-"NA"
      }else{
        y=t.test(as.numeric(a[i,(1:rep)]),as.numeric(a[i,(rep+1):(rep*2)]))
        Pvalue[i]<-y$p.value
        log2_FC[i]<-log2((mean(as.numeric(a[i,(1:rep)]))+0.001)/(mean(as.numeric(a[i,(rep+1):(rep*2)]))+0.001)) 
        fdr[i]=p.adjust(Pvalue[i], "BH") 
      }
    }
    # 在原文件后面加入log2FC，p value和FDR,共3列；
    out<-cbind(a,log2_FC,Pvalue,fdr)
    
    # out$tax=otu_table$compound
    head(out)
    x = out
    
    WT <-subset(out,fdr < 0.05 & log2_FC != "NA")
    dim(WT)
    head(WT)
    
    ###########添加物种丰度
    # dim(count)
    # str(count)
    count = as.matrix(count)
    norm = t(t(count)/colSums(count)) #* 100 # normalization to total 100
    dim(norm)
    norm1 = norm %>% 
      t() %>% as.data.frame()
    # head(norm1)
    #数据分组计算平均值
    library("tidyverse")
    dim(norm1)
    
    iris.split <- split(norm1,as.factor(sub_design$Group))
    iris.apply <- lapply(iris.split,function(x)colMeans(x))
    # 组合结果
    iris.combine <- do.call(rbind,iris.apply)
    norm2= t(iris.combine)
    
    #head(norm)
    str(norm2)
    norm2 = as.data.frame(norm2)
    # dim(x)
    # dim(norm2)
    x = cbind(x,norm2)
    # head(x)
    #在加入这个文件taxonomy时，去除后面两列不相干的列
    # 读取taxonomy，并添加各列名称
    vegan_tax <-  function(physeq){
      tax <-  tax_table(physeq)
      
      return(as(tax,"matrix"))
    }
    taxonomy = as.data.frame(vegan_tax(ps_rela))
    head(taxonomy)
    # taxonomy <- as.data.frame(tax_table(ps1))
    
    #发现这个注释文件并不适用于直接作图。
    #采用excel将其分列处理，并且删去最后一列，才可以运行
    if (length(colnames(taxonomy)) == 6) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    }else if (length(colnames(taxonomy)) == 7) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
    }else if (length(colnames(taxonomy)) == 8) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","rep")
    }
    # colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    
    # Taxonomy排序，并筛选OTU表中存在的
    library(dplyr)
    taxonomy$id=rownames(taxonomy)
    # head(taxonomy)
    tax = taxonomy[row.names(x),]
    
    if (length(colnames(taxonomy)) == 7) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      # x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 8) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 9) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }
    
    head(x)
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    
    filename = paste(path,"/a3_Desep_for_",group,"_allOTU_.csv",sep = "")
    write.csv(x,filename,quote = FALSE)
    
    filename = paste(path,"/a3_Desep_for_",group,"differ_OTU.csv",sep = "")
    write.csv(WT,filename,quote = FALSE)
    
    
  }
  
}


# Wlx(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
Wlx = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv"){
  library(phyloseq)
  #导入otu表格
  otu = read.delim(otu,row.names = 1)
  head(otu)
  otu = as.matrix(otu)
  str(otu)
  #导入注释文件
  tax = read.delim(tax,row.names = 1)
  head(tax)
  tax = as.matrix(tax)
  # taxa_names(tax)
  
  #导入分组文件
  map = read.delim(map,row.names = 1)
  head(map)
  
  # #导入进化树
  # tree = read.tree("./otus.tree")
  # tree
  
  ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                 sample_data(map) ,
                 tax_table(tax)
                 # phy_tree(tree)
                 
  )
  ps
  
  ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  sub_design <- as.data.frame(sample_data(ps))
  
  sub_design <- as.data.frame(sample_data(ps))
  Desep_group <- as.character(levels(sub_design$Group))
  Desep_group
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  otu_table = as.data.frame(vegan_otu(ps))
  
  count = as.matrix(otu_table)
  count <- t(count)
  #数据整理形式一######
  
  sub_design <- as.data.frame(sample_data(ps))
  sub_design $ID = row.names(sub_design )
  dim(sub_design)
  sub_design$SampleType = as.character(sub_design$Group)
  sub_design$SampleType <- as.factor(sub_design$SampleType)
  
  aaa = combn(Desep_group,2)
  ##设置文件路径
  path = "./Wilcox"
  
  dir.create(path)
  # aaa[,1]
  i = 1
  for (i in 1:dim(aaa)[2]) {
    Desep_group = aaa[,i]
    # 转换原始数据为百分比，
    norm = t(t(count)/colSums(count,na=TRUE)) * 100 # normalization to total 100
    head(norm)
    norm=as.data.frame(norm)
    a=norm
    head(a)
    #预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
    Pvalue<-c(rep(0,nrow(a)))
    fdr<-c(rep(0,nrow(a)))
    log2_FC<-c(rep(0,nrow(a)))
    
    library(tidyverse)
    df_filter<- filter(sub_design ,SampleType %in% Desep_group)
    
    head(df_filter)
    
    a = as.data.frame(a)
    a = a[as.character(df_filter$ID)]
    head(a)
    rep = length(as.character(df_filter$ID))/2
    
    a = as.matrix(a)
    ###########开始运行脚本
    for(i in 1:nrow(a)){
      if(sd(a[i,(1:rep)])==0&&sd(a[i,(rep+1):(rep*2)])==0){
        Pvalue[i] <-"NA"
        log2_FC[i]<-"NA"
      }else{
        y=t.test(as.numeric(a[i,(1:rep)]),as.numeric(a[i,(rep+1):(rep*2)]))
        Pvalue[i]<-y$p.value
        log2_FC[i]<-log2((mean(as.numeric(a[i,(1:rep)]))+0.001)/(mean(as.numeric(a[i,(rep+1):(rep*2)]))+0.001)) 
        fdr[i]=p.adjust(Pvalue[i], "BH") 
      }
    }
    
    for(i in 1:nrow(a)){
      if(sd(a[i,(1:rep)])==0&&sd(a[i,(rep+1):(rep*2)])==0){
        Pvalue[i] <-"NA"
        log2_FC[i]<-"NA"
      }else{
        y=wilcox.test(as.numeric(a[i,(1:rep)]),as.numeric(a[i,(rep+1):(rep*2)]),exact=FALSE)
        Pvalue[i]<-y$p.value
        log2_FC[i]<-log2((mean(as.numeric(a[i,(1:rep)]))+0.001)/(mean(as.numeric(a[i,(rep+1):(rep*2)]))+0.001)) 
        fdr[i]=p.adjust(Pvalue[i], "BH") 
      }
    }
    
    # 在原文件后面加入log2FC，p value和FDR,共3列；
    out<-cbind(a,log2_FC,Pvalue,fdr)
    
    # out$tax=otu_table$compound
    head(out)
    x = out
    
    WT <-subset(out,fdr < 0.05 & log2_FC != "NA")
    dim(WT)
    head(WT)
    
    ###########添加物种丰度
    # dim(count)
    # str(count)
    count = as.matrix(count)
    norm = t(t(count)/colSums(count)) #* 100 # normalization to total 100
    dim(norm)
    norm1 = norm %>% 
      t() %>% as.data.frame()
    # head(norm1)
    #数据分组计算平均值
    library("tidyverse")
    dim(norm1)
    
    iris.split <- split(norm1,as.factor(sub_design$SampleType))
    iris.apply <- lapply(iris.split,function(x)colMeans(x))
    # 组合结果
    iris.combine <- do.call(rbind,iris.apply)
    norm2= t(iris.combine)
    
    #head(norm)
    str(norm2)
    norm2 = as.data.frame(norm2)
    # dim(x)
    # dim(norm2)
    x = cbind(x,norm2)
    # head(x)
    #在加入这个文件taxonomy时，去除后面两列不相干的列
    # 读取taxonomy，并添加各列名称
    vegan_tax <-  function(physeq){
      tax <-  tax_table(physeq)
      
      return(as(tax,"matrix"))
    }
    taxonomy = as.data.frame(vegan_tax(ps))
    head(taxonomy)
    # taxonomy <- as.data.frame(tax_table(ps1))
    
    #发现这个注释文件并不适用于直接作图。
    #采用excel将其分列处理，并且删去最后一列，才可以运行
    if (length(colnames(taxonomy)) == 6) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    }else if (length(colnames(taxonomy)) == 7) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
    }else if (length(colnames(taxonomy)) == 8) {
      colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","rep")
    }
    # colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus")
    
    # Taxonomy排序，并筛选OTU表中存在的
    library(dplyr)
    taxonomy$id=rownames(taxonomy)
    # head(taxonomy)
    tax = taxonomy[row.names(x),]
    
    if (length(colnames(taxonomy)) == 7) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      # x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 8) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }else if (length(colnames(taxonomy)) == 9) {
      x = x[rownames(tax), ] # reorder according to tax
      x$phylum = gsub("","",tax$phylum,perl=TRUE) 
      x$class = gsub("","",tax$class,perl=TRUE) 
      x$order = gsub("","",tax$order,perl=TRUE) 
      x$family = gsub("","",tax$family,perl=TRUE) 
      x$genus = gsub("","",tax$genus,perl=TRUE) 
      x$species = gsub("","",tax$species,perl=TRUE) 
    }
    
    head(x)
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    
    filename = paste(path,"/a3_Desep_for_",group,"_allOTU_.csv",sep = "")
    write.csv(x,filename,quote = FALSE)
    
    filename = paste(path,"/a3_Desep_for_",group,"differ_OTU.csv",sep = "")
    write.csv(WT,filename,quote = FALSE)
    
    
  }
  
}


# Diff_test (otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",method ="DESep2")
Diff_test = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",method ="DESep2" ){
  library(phyloseq)
  if (method == "DESep2") {
    DESep2(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
  }
  if (method == "Edger") {
    Edger(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
    
  }
  if (method == "Ttest") {
    Ttest(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
    
  }
  if (method == "Wlx") {
    Wlx(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
    
  }
  
}
