# 清空工作空间
rm(list=ls())
# 安装R包
install.packages("seqinr")
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

# 第一步：导出抬头格式 
geshi <- c("#mega","!Title:;")
write.table(geshi,"dist.meg",quote = F,row.names = F,col.names = F)

# 第二步：导出序列描述
r_name <- paste0("[",1:35,"]")

description <- seqinr::getAnnot(S)
description <- gsub(">","#",description) # 将>批量转换为#
description <- gsub(" ","_",description) # 将空格批量转换为_
description <- data.frame(description)
row.names(description) <- r_name
write.table(description,"dist.meg",append = T,quote = F,col.names = F)

# 第三步：导出矩阵的列抬头
c_name <- c("[",1:35,"]")
write.table(t(c_name),"dist.meg",append = T,quote = F,row.names = F,col.names = F,sep = "\t")

row.names(tz) <- r_name # 将行名转为[]
# 用欧氏距离定义样本间的距离
dt <- dist(scale(tz))
dt <- data.matrix(dt)
dt[upper.tri(dt,diag = T)] <- NA # 上三角矩阵和主对角线都赋值为NA

# 第四步：导出矩阵
write.table(round(dt,4),"dist.meg",append = T,quote = F,col.names = F,na="",sep = "\t")
# 将矩阵写入到文件中，让NA值显示为空，间隔为制表符
