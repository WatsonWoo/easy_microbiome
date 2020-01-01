
# #清空内存
# rm(list=ls()) 
# # network
# otu = "./otutab.txt"
# tax = "./taxonomy.txt"
# map = "./metadata.tsv"
# N = 0.05
# 
# r.threshold=0.6
# p.threshold=0.05
# label = TRUE


#network (otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",N = 0.05,r.threshold=0.6,p.threshold=0.05,label = TRUE)
network = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",N = 0.05,r.threshold=0.6,p.threshold=0.05,label = FALSE){
  #注意使用原始coun丰度文件
  path = "./network/"
  dir.create(path)
  library("phyloseq")
  library(ggcorrplot)
  library(igraph)
  library(psych)
  library(network)
  library(ggplot2)
  library(sna)
  library(ergm)
  library(igraph)
  
  
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
  
  
  
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)
    
    return(as(tax,"matrix"))
  }
  
  
  #相对丰度转换
  ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) )
  ps_sub = filter_taxa(ps_rela, function(x) sum(x ) > N , TRUE)#筛选序列数量大于1的
  ps_sub
  
  
  
  
  #读入分组文件
  design = mapping= as.data.frame(sample_data(ps_sub))
  head(design)
  
  ### 提取丰度表格
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  head(otu_table)
  
  
  
  
  # 提取注释文件表格
  tax_table = as.data.frame(vegan_tax(ps_sub ))
  head(tax_table)
  #合并丰度和
  
  # library("corrplot")
  
  # corr.test求相关性矩阵
  
  network = otu_table
  dim(network)
  
  y = matrix(1:2409,nrow = 14,ncol = length(unique(design$Group)))
  d = N
  
  
  layouts = as.character(unique(design$Group))
  library("tibble")
  head(t(network))
  
  mapping$ID = row.names(mapping)
  ##################---------------------------------------开始计算网络-----------------------------------------------------------------
  
  # 加载网络性质的计算函数
  source("./net_pro.R")
  # 加载节点性质的计算函数
  source("./node_pro.R")
  ########相关函数加载#########
  aa = 1
  plots = list()
  layout = layouts[1]
  
  for (layout in layouts) {
    # network_sub = iris.split[[layout]]
    
    
    library(tidyverse)
    
    xx = dplyr::filter(mapping, Group %in% layout)
    network_sub = network[xx$ID]
    head(network_sub)
    
    dim(network_sub)
    n=ncol(network_sub)
    #
    network_sub[n+1]=apply(network_sub[c(1:nrow(network_sub)),],1,sum)
    #
    network_sub=network_sub[network_sub[n+1] > 0,1:n]
    
    head(network_sub)
    
    
    occor = corr.test(t(network_sub),use="pairwise",method="pearson",adjust="fdr",alpha=.05)
    # occor = cor(t(network_sub),use="pairwise",method="spearman")
    
    occor.r = occor$r # 取相关性矩阵R值
    
    occor.p = occor$p # 取相关性矩阵p值
    
    occor.r[occor.p > p.threshold|abs(occor.r)<r.threshold] = 0
    
    head(occor.r)
    
    dim(occor.r)
    g <- network::network(occor.r, directed=FALSE)
    (summary(g))
    g$gal
    
    net  = g
    m <- as.matrix.network.adjacency(net)  # get sociomatrix
    
    
    
    
    # get coordinates from Fruchterman and Reingold's force-directed placement
    # algorithm.
    plotcord <- data.frame(gplot.layout.fruchtermanreingold(m, NULL))
    colnames(plotcord) = c("X1", "X2")
    dim(plotcord)
    plotcord$elements <- colnames(occor.r)
    edglist <- as.matrix.network.edgelist(net)
    edglist = as.data.frame(edglist)
    head(edglist)
    #
    set.edge.value(g,"weigt",occor.r)
    net  = g
    # ?get.edge.attribute
    
    
    edglist$weight = as.numeric(network::get.edge.attribute(net,"weigt"))
    dim(edglist)
    edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
    head(edges)
    edges$weight = as.numeric(network::get.edge.attribute(net,"weigt"))
    ##这里将边权重根据正负分为两类
    aaa = rep("a",length(edges$weight))
    for (i in 1:length(edges$weight)) {
      if (edges$weight[i]> 0) {
        aaa[i] = "+"
      }
      if (edges$weight[i]< 0) {
        aaa[i] = "-"
      }
    }
    #添加到edges中
    edges$wei_label = as.factor(aaa)
    colnames(edges) <- c("X1", "Y1","OTU_1", "X2", "Y2","OTU_2","weight","wei_label")
    edges$midX <- (edges$X1 + edges$X2)/2
    edges$midY <- (edges$Y1 + edges$Y2)/2
    head(edges)
    
    dim(plotcord)
    ##plotcord这个表格进行物种注释的添加
    row.names(plotcord) = plotcord$elements
    dim(plotcord)
    
    head(plotcord)
    
    library(ggrepel)
    
    network_tax =tax_table
    head(network_tax)
    
    
    
    dim(network_tax)
    
    
    
    res = merge(plotcord,network_tax,by = "row.names",all = F)
    dim(res)
    head(res)
    row.names(res) = res$Row.names
    res$Row.names = NULL
    plotcord = res
    head(plotcord)
    plotcord$mean  =rowMeans(otu_table)
    
    pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                    data = edges, size = 0.5) +
      geom_point(aes(X1, X2,size = mean,fill = "Phylum" ),pch = 21, data = plotcord) + scale_colour_brewer(palette = "Set1") +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      labs( title = paste(layout,"network",sep = "_"))+
      # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
      # discard default grid + titles in ggplot2
      theme(panel.background = element_blank()) +
      # theme(legend.position = "none") +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",  colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    
    if (label == TRUE ) {
      pnet <- pnet +  geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)
    }
    pnet
    plotname = paste(path,"/network",layout,".pdf",sep = "")
    
    
    ggsave(plotname, pnet, width = 12, height =8)
    
    
    plots[[aa]] = pnet
    
    occor.r = as.matrix(occor.r)
    # # 构建igraph对象构建邻接矩阵
    igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
    igraph
    ###网络边的赋值及其设置
    igraph.weight <- E(igraph)$weight# 将igraph weight属性赋值到igraph.weight,用于后边做图
    E(igraph)$weight <- NA
    igraph<-remove.edge.attribute(igraph,"weight")#把边值删除
    netpro_result<-net_pro(igraph)
    colnames(netpro_result)<-layout
    y = as.data.frame(y)
    colnames(y) = layouts
    # head(y)
    y[layout] = netpro_result[,1]
    row.names(y) = row.names(netpro_result)
    y
    
    aa = aa+1
  }
  
  
  
  plotname = paste(path,"/network_all.pdf",sep = "")
  
  library(ggpubr)
  p  = ggarrange(plotlist = plots,ncol=4, nrow=2, common.legend = TRUE, legend="right")
  # p
  
  ggsave(plotname, p,width = 18,height = 10)
  
  plotname1 = paste(path,"/network_all.jpg",sep = "")
  ggsave(plotname1, p,width = 18,height = 10)
  
  
  tablename <- paste(path,"/co-occurrence_net_net_pro4",".csv",sep = "")
  write.csv(y,tablename)
  
  
  #下面多中可视化方法共选择
  # plotcord <- data.frame(gplot.layout.kamadakawai(m, NULL))
  # plotcord <- data.frame(gplot.layout.adj(m, NULL))
  # plotcord <- data.frame(gplot.layout.circle(m, NULL))
  # plotcord <- data.frame(gplot.layout.circrand(m, NULL))
  # plotcord <- data.frame(gplot.layout.eigen(m, NULL))
  # plotcord <- data.frame(gplot.layout.geodist(m, NULL))
  # plotcord <- data.frame(gplot.layout.hall(m, NULL))
  # plotcord <- data.frame(gplot.layout.mds(m, NULL))
  # plotcord <- data.frame(gplot.layout.princoord(m, NULL))
  # plotcord <- data.frame(gplot.layout.random(m, NULL))
  # plotcord <- data.frame(gplot.layout.rmds(m, NULL))
  # plotcord <- data.frame(gplot.layout.segeo(m, NULL))
  # plotcord <- data.frame(gplot.layout.seham(m, NULL))
  # plotcord <- data.frame(gplot.layout.spring(m, NULL))
  # plotcord <- data.frame(gplot.layout.springrepulse(m, NULL))
  # plotcord <- data.frame(gplot.layout.target(m, NULL))
  
  
  
}


