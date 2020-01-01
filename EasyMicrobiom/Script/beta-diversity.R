# beta-diversity

#beta(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv")
beta <- function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv"){
    path = "./bata/"
    dir.create(path)
    library(phyloseq)
    library(vegan)
    library(ggplot2)
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
    
    ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela 
    unif <- distance(ps1_rela , method="bray", type="samples")
    #这里请记住pcoa函数
    pcoa = cmdscale(unif, k=2, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
    
    points = as.data.frame(pcoa$points) # 获得坐标点get coordinate string, format to dataframme
    colnames(points) = c("x", "y") #命名行名
    eig = pcoa$eig
    #eig特征值得到
    sub_design = as.data.frame(sample_data(ps1_rela))
    points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
    #write.table(points,"pcoa_bray_curtis.txt",quote = FALSE,row.names = F,
    #           col.names = T,sep = "\t")
    head(points)
    
    ado = adonis(unif~ sub_design$Group,permutations = 999,method="bray")
    a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 <- paste("adonis:R ",a, sep = "")
    b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",b, sep = "")
    title = paste(R2," ",p_v, sep = "")
    title
    
    # mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A","#E6AB02", "#B3DE69")
    # jpeg(file="./a2_bray_PCOA.jpeg")
    p2 <-ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
           title=title)+
      stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))+
      #stat_ellipse( linetype = 1,level = 0.8)+
      #geom_text_repel(aes(label=points$id),size=4)+
      # scale_colour_manual(values = mi,guide = guide_legend(title = NULL))+
      # scale_fill_manual(values = mi,guide = guide_legend(title = NULL))+
      #labs(title = "toamto hea and dis")+
      guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL)) 
    # p2
    # points$id=row.names(points)
    # p+geom_text(aes(label=points$id),size=4)#?stat_ellipse
    p2 = p2+theme_bw()+
      
      #scale_y_continuous(expand = c(0,0))+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
      theme(
        
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
        #legend.position = "none"#是否删除图例
        
      ) 
    p2
    FileName <- paste(path,"./a2_bray_PCOA.pdf", sep = "")
    ggsave(FileName, p2, width = 12, height = 8)
    FileName1 <- paste(path,"./a2_bray_PCOA.jpg", sep = "")
    ggsave(FileName1 , p2, width = 12, height = 8)
    
  
    
  }
