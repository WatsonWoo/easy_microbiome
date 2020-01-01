# 计算alpha多样性和出图
#alpha(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",inde="Shannon",Plot = TRUE )

alpha = function(otu,tax,map,inde="Shannon",Plot = TRUE){
  path = "./alpha/"
  dir.create(path)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(picante)      #picante 包加载时默认同时加载 vegan
  library(agricolae)
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
  
  
  
  ##按照最小序列数抽平
  total = mean(sample_sums(ps));total
  # total = min(sample_sums(ps));total
  standf = function(x,t = total)round(t*(x/sum(x)))
  ps11 = transform_sample_counts(ps,standf)
  
  mapping = sample_data(ps11)
  
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  count = as.data.frame(t(vegan_otu(ps11)))
  head(count)
  
  # 加载VEGAN包
  alpha=diversity(count, "shannon")
  
  x = t(count) ##转置，行为样本，列为OTU
  head(x)
  
  
  ##核心，计算两个指标
  Shannon = vegan::diversity(x)  ##默认为shannon
  Shannon
  Inv_Simpson <- vegan::diversity(x, index = "invsimpson")
  Inv_Simpson
  
  #计算OTU数量
  S <- vegan::specnumber(x);S  ##每个样本物种数。等价于S2 = rowSums(x>0)
  S2 = rowSums(x>0)
  
  
  #多样性指标：均匀度Pielou_evenness，Simpson_evenness
  Pielou_evenness <- Shannon/log(S)
  Simpson_evenness <- Inv_Simpson/S
  
  
  
  
  
  
  
  est <- estimateR(x)
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  
  report = cbind(Shannon, Inv_Simpson, Pielou_evenness, Simpson_evenness,
                 Richness, Chao1,ACE) #不同列进行合并
  head(report)
  
  index = merge(mapping,report , by="row.names",all=F)
  head(index)
  FileName <- paste(path,"./alpha_diversity.csv", sep = "")
  write.csv(index,FileName,sep = "")
  
  # 
  # if (plot == FALSE ) {NULL}
  if (Plot == TRUE ) {
    Mytheme <- theme_bw()+
      
      # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
      theme(
        
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold"),
        legend.position = "none"#是否删除图例
        
      )
    
    # head(data_wt )
    data_wt = index
    data_wt$group = data_wt$Group 
    # alpha = "Shannon"
    # data_wt = as.matrix(data_wt)
    data_wt[inde]
    
    #构造待分析的子数据框
    ss <- data_wt[inde]
    colnames(ss) <- c("count")
    ss$group = data_wt$group
    
    
    name_i = colnames(data_wt[inde])
    #求取均值和方差
    wen1 = as.data.frame(tapply(as.vector(as.matrix(data_wt[inde])),data_wt$group,mean,na.rm=TRUE))
    wen2 = as.data.frame(tapply(as.vector(as.matrix(data_wt[inde])),data_wt$group,sd,na.rm=TRUE))
    went = cbind(wen1,wen2)
    
    colnames(went) = c("mean" ,"SD")
    went
    
    #进行方差检验 下面wtx3为提取的p值
    model<-aov(count ~ group, data= ss)#方差分析
    wtx1 = summary(model)
    wtx2 = wtx1[[1]]
    wtx3 = wtx2[5]#
    
    
  
  
    # out <- LSD.test(model,"group", p.adj="none")#进行多重比较，不矫正P值
    # aa = out$group#结果显示：标记字母法
    # aa$group = row.names(aa)
    # aa
    # ? duncan.test
    out <-duncan.test (model,"group")
    aa = out$groups# 查看每个组的label
    
    aa$group = row.names(aa)
    stat = aa
    aa
    
    
    wentao = merge(aa,went, by="row.names",all=F)
    wentao
    # FileName <- paste(name_i,"_aov_bar", ".csv", sep = "_")
    # write.csv(wentao,FileName,quote = F)
    # colnames(wentao) = c(colnames(wentao[1:4]),"mean" ,"SD")
    #使用的tidyverse函数，对数据框添加两列，目的为柱状图添加bar
    aa = mutate(wentao, ymin = mean - SD, ymax =  mean + SD)
    a = max(aa$mean)*1.5##用于设置y轴最大值
    
    ### 出图柱状图
    p = ggplot(aa , aes(x = group, y = mean,colour= group)) +
      geom_bar(aes(colour= group,fill = group),stat = "identity", width = 0.4,position = "dodge") +
      
      geom_errorbar(aes(ymin=ymin,
                        ymax=ymax),
                    colour="black",width=0.1,size = 1)+
      # geom_hline(aes(yintercept=mean(as.vector(as.matrix(data_wt[i])))), colour="black", linetype=2) +
      # geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      scale_y_continuous(expand = c(0,0),limits = c(0,a))+
      labs(x=paste(name_i,"", sep = ""),
           y =name_i )
    p
    
    
    p = p + geom_text(aes(label = groups,y=ymax, x = group,vjust = -0.3,size = 6))
    p
    
    #as.vector(as.matrix(data_wt[i]))为进行差异分析的一组数据
    p=p+Mytheme
    p
    
    if (length(unique(data_wt$group))>3){	p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
    p
    FileName <- paste(path,name_i,"_aov_bar", ".pdf", sep = "")
    # library("Cairo")
    ggsave(FileName, p, width = 10, height = 8)
    FileName1 <- paste(path,name_i,"_aov_bar", ".jpg", sep = "")
    # library("Cairo")
    ggsave(FileName1, p, width = 10, height = 8)
    
  }
}






