# 清空工作空间
rm(list=ls())
# 安装R包
install.packages("ggplot2")
install.packages("factoextra")
install.packages("seqinr")
install.packages("dplyr")

# 加载13维氨基酸理化指标
index <- read.table("physicochemical_properties.txt")
# 将行名全部修改成大写
row.names(index) <- toupper(row.names(index))

# 计算距离矩阵
dt <- dist(scale(index[,c(8,7,10,11)]))

#===================系统聚类===================
# 系统聚类
hc <- hclust(dt,method = "average")
# 展示聚类谱系图
library(ggplot2);library(factoextra)
graph0 <- fviz_dend(hc, 
            k=6,                 ## 分成6类
            cex = 1,              ## 设置数据标签的字体大小
            horiz = T,                ## 水平摆放图形
            k_colors = 2:8,
            labels_track_height=0,
            color_labels_by_k=T,    ## 自动设置数据标签颜色
            lwd = 1               ## 设置分支的线宽
            ) +      ## 绘制矩形聚类图
  theme(axis.title.x = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.8)))  ## 修改坐标字体大小
graph0
# 加载序列数据
S <- seqinr::read.fasta(file.choose())
# 提取特征向量
get_tz <- function(m){
  lst <- list()
  for (i in 1:m) {
    n <- numeric(6)
    weizhi <- numeric(6)
    for (j in 1:length(S[[i]])) {
      s <- S[[i]][j]
      if (s %in% c("a","i","l","m","v","f")){
        n[1]=n[1]+1
        weizhi[1]=weizhi[1]+j
      }
      if (s %in% c("y")){
        n[2]=n[2]+1
        weizhi[2]=weizhi[2]+j
      }
      if (s %in% c("q","n","s","g","t","p","w")){
        n[3]=n[3]+1
        weizhi[3]=weizhi[3]+j
      }
      if (s %in% c("k","r","h")){
        n[4]=n[4]+1
        weizhi[4]=weizhi[4]+j
      }
      if (s %in% c("e","d")){
        n[5]=n[5]+1
        weizhi[5]=weizhi[5]+j
      }
      if (s %in% c("c")){
        n[6]=n[6]+1
        weizhi[6]=weizhi[6]+j
      }
    }
    
    lst[[i]] <- n/length(S[[i]])
    lst[[i]][7:12] <- weizhi/n
  }
  return(lst)
}

# 运行函数所需时间
system.time(tz <- get_tz(length(S)))
tz <- get_tz(length(S)) |> data.frame() |> t()

library(dplyr)
description <- seqinr::getAnnot(S) %>%  
  gsub(">","",.) # 将>批量转换为空格

row.names(tz) <- description
#===================相似性距离矩阵===================
# 计算距离矩阵
dt <- dist(scale(tz))
dt_show <- as.matrix(round(dt,3))
#===================系统聚类===================
# 系统聚类
hc <- hclust(dt,method = "average")

# 展示聚类谱系图
graph <- fviz_dend(hc, 
                   k=5,                 ## 分成5类
                   cex = 0.7,              ## 设置数据标签的字体大小
                   horiz=T,                ## 水平摆放图形
                   k_colors = c(2:4,6,7),
                   labels_track_height=5,
                   ylab = "Linkage distance",
                   color_labels_by_k=T,    ## 自动设置数据标签颜色
                   lwd = 1,               ## 设置分支的线宽
                   type = "rectangle") +      ## 绘制矩形聚类图
  theme(axis.title.y = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1))) + ## 修改坐标字体大小
  labs(y="") ## 环形进化树不用显示y

graph

