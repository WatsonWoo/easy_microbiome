# otu = "./otutab.txt"
# tax = "./taxonomy.txt"
# map = "./metadata.tsv"
# rep = 6
# #清空内存
# rm(list=ls()) 




#flowerplot(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",rep = 6)

flowerplot = function(otu = "./otutab.txt",tax = "./taxonomy.txt",map = "./metadata.tsv",rep = 6){
  
  path = "./flower_Ven/"
  dir.create(path)
  library("phyloseq")
  library(UpSetR)
  library("tibble")
  library (VennDiagram) 
  
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
  mapping = as.data.frame(sample_data(ps))
  ps1 = ps
  ps1
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  aa = vegan_otu(ps1)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  dim(count)
  sub_design <- as.data.frame(sample_data(ps1))
  
  # levels(sub_design$SampleType)[1]
  # name1 = paste("name",levels(sub_design$SampleType)[1],sep = "")
  sub_design $SampleType = gsub("1","",sub_design $Group)
  sample_data(ps ) = sub_design
  levels(sub_design $Group)
  
  ##########这里的操作为提取三个分组
  pick_val_num <- rep *2/3
  count[count > 0] <- 1###这个函数只能用于0,1 的数据，所以我这么转换
  
  count2 = as.data.frame(count)
  
  library("tibble")
  #数据分组
  iris.split <- split(count2,as.factor(sub_design$Group))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  head(ven2)
  ven2[ven2 < pick_val_num]  = 0
  ven2[ven2 >=pick_val_num]  = 1
  ven2 = as.data.frame(ven2)
  
  
  #########更加高级的设置在这里可以查看#https://mp.weixin.qq.com/s/6l7gftKQfiyxNH66i19YtA
  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)
  head(ven2)
  
  
  all_num = dim(ven2[rowSums(ven2) == 4,])[1]
  
  ven2[,1] == 1
  
  A = rep("A",length(colnames(ven2)))
  B = rep(1,length(colnames(ven2)))
  
  i = 1
  for (i in 1:length(colnames(ven2))) {
    B[i] = length(ven2[rowSums(ven2) == 1,][,i][ven2[rowSums(ven2) == 1,][,i] == 1])
    
    A[i] = colnames(ven2)[i]
  }
  A
  B
  
  library("plotrix")
  
  flower_plot <- function(sample, value, start, a, b, cer_label,col_c,
                          ellipse_col = rgb(135, 206, 235, 150, max = 255),
                          circle_col = rgb(0, 162, 214, max = 255),
                          circle_text_cex = 1.5
  ) {
    par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
    plot(c(0,10),c(0,10),type="n")
    n   <- length(sample)
    deg <- 360 / n
    res <- lapply(1:n, function(t){
      draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                   y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                   col = ellipse_col,
                   border = ellipse_col,
                   a = a, b = b, angle = deg * (t - 1))
      text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
           value[t]
      )
      
      if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) - start,
             adj = 1,
             cex = circle_text_cex
        )
        
      } else {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) + start,
             adj = 0,
             cex = circle_text_cex
        )
      }			
    })
    draw.circle(x = 5, y = 5, r = 1.3, col = circle_col, border = circle_col)
    
    text(x = 5, y = 5,cer_label,cex = 2,col = col_c)
  }
  
  
  library(RColorBrewer)#调色板调用包
  
  # #调用所有这个包中的调色板
  # display.brewer.all()
  # #提取特定个数的调色板颜色，会出图显示
  # display.brewer.pal(9,"Oranges")
  # display.brewer.pal(8,"Set1")
  # ##仅仅只显示色号,我们要在颜色上并显示色号才好
  mi = brewer.pal(9,"Blues")
  mi = brewer.pal(9,"Set1")
  
  # library("scales")
  # show_col(mi)
  
  
  
  
  pdf(paste(path,"/flower.pdf",sep = ""), width = 6, height = 6)
  
  flower_plot(A,
              B,
              ellipse_col =  mi[6],
              circle_col = mi[4],
              90, 0.4, 2,
              cer_label = paste("all",":",all_num,sep = ""), col_c = mi[1]
  )
  
  
  
  dev.off()
  
  jpeg(paste(path,"/flower.jpg",sep = ""))
  
  flower_plot(A,
              B,
              ellipse_col =  mi[6],
              circle_col = mi[4],
              90, 0.4, 2,
              cer_label = paste("all",":",all_num,sep = ""), col_c = mi[1]
  )
  
  
  
  dev.off()
}



