
# otu = "./otutab.txt"
# tax = "./taxonomy.txt"
# map = "./metadata.tsv"

# 接下来我门绘制门水品的物种词云
# j = "Phylum"
# j = "Class"
# j = "Order"
# j =  "Family"
# j = "Genus"

# 微生物生态-扩增子数据词云



# 基础词云，用于描述不同分组就某个分类水平绘制词云


result = wordTax(otu = "./otutab.txt",tax = "./taxonomy.txt", map = "./metadata.tsv",j = "Class")

# 提取图片
p1 = result [[1]]
p1
# 提取作图数据
data = result [[2]]

#
path = getwd()
FileName2 <- paste(path,"/word—base.pdf", sep = "")

ggsave(FileName2, p1, width = 8, height =6)



# facet模式，用于按照门水平绘制词云，也就是是纵坐标按照门水平分类，横坐标按照分组展示。

result = wordTax(otu = "./otutab.txt",tax = "./taxonomy.txt", map = "./metadata.tsv",j = "Family",facet = TRUE)

# 提取图片
p1 = result [[1]]

# 提取作图数据
data = result [[2]]

#
path = getwd()
FileName2 <- paste(path,"/word—facet.pdf", sep = "")

ggsave(FileName2, p1, width = 12, height =50,limitsize = FALSE)












wordTax = function(otu = "./otutab.txt",tax = "./taxonomy.txt", map = "./metadata.tsv",j = "Phylum",facet = FALSE){
  
  
  
  library(tidyverse)
  library(ggwordcloud)
  library("phyloseq")
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
  
  
  
  
  ps_sub <- ps %>% tax_glom(taxrank = j) %>% transform_sample_counts(function(x) {x/sum(x)} )
  ps_sub 
  otu_table = as.data.frame(t(vegan_otu(ps_sub  )))
  
  head(otu_table)
  
  # 按照不同分组求取每组平均丰度列表
  count = t(otu_table)
  count2 = as.data.frame(count)
  head(count)
  #提取分组文件
  sub_design <- as.data.frame(sample_data(ps))
  
  #构造分组
  aa = sub_design[,"Group"]
  colnames(aa) = "Vengroup"
  iris.split <- split(count2,as.factor(aa$Vengroup))
  
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colMeans(x[]))
  
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  # sum(ven2[,1])#验证一下是否正确
  ven2 = as.data.frame(ven2)
  head(ven2)
  ven2$ID  =row.names(ven2)
  # 提取物种注释文件：
  
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)
    
    return(as(tax,"matrix"))
  }
  tax_table = as.data.frame(vegan_tax(ps_sub))
  head(tax_table)
  
  # match(row.names(tax_table),ven2$ID) 这是匹配的，所以直接合并
  
  ven3 = cbind(tax_table,ven2)
  
  
  #数据宽边长
  library(reshape2)
  plotdata<-melt(ven3,id.vars = c("ID",colnames(tax_table)),
                 variable.name='Group',
                 value.name='mean')
  head(plotdata)
  
  set.seed(42)
  p = ggplot(plotdata, aes_string(label = j, size = "mean",color = "Phylum")) +
    geom_text_wordcloud() +
    scale_size_area(max_size = 15) +
    theme_minimal() +
    facet_wrap(~Group)
  p
  
  if(facet == TRUE){
    
    p = ggplot(plotdata, aes_string(label = j, size = "mean",color = "Phylum")) +
      geom_text_wordcloud() +
      scale_size_area(max_size = 15) +
      theme_minimal() +facet_grid(Phylum~Group)
    # facet_wrap(~Group)
    p
  }
  

  
  return(list(p,plotdata))
}

# 添加数量标记
# ggplot(plotdata, aes_string(label = j, size = "mean",color = "Phylum",
#                             label_content = sprintf("%s<span style='font-size:7.5pt'>(%g)</span>")  
#                             )) +
#   geom_text_wordcloud() +
#   scale_size_area(max_size = 15) +
#   theme_minimal() 




















# #设置单词角度
# angle = 45 * sample(-2:2, nrow(plotdata),replace = TRUE,prob = c(1, 1, 4, 1, 1))
# ggplot(plotdata, aes(label = ID, size = mean,color = Phylum)) +
#   geom_text_wordcloud(aes(angle = angle )) +
#   scale_size_area(max_size = 15) +
#   theme_minimal() +
#   facet_wrap(~Group)

# set.seed(42)
# ggplot(plotdata, aes(label = ID, size = mean,color = Phylum)) +
#   geom_text_wordcloud(aes(angle = 45 * sample(-2:2, nrow(plotdata),
#                                               replace = TRUE,prob = c(1, 1, 4, 1, 1)))) +
#   scale_size_area(max_size = 15) +
#   theme_minimal() +
#   facet_wrap(~Group)
#   










