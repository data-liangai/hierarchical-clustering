# 清空工作空间
import gc
gc.collect()

# 安装需要的Python库
# pip install pandas numpy scipy matplotlib seaborn biopython

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import zscore
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

os.chdir("D:\hierarchical-clustering")
# 加载13维氨基酸理化指标
index = pd.read_csv("physicochemical_properties.txt", sep="\t", index_col=0)
# 将行名全部修改成大写
index.index = index.index.str.upper()
# 计算距离矩阵之前进行标准化
scaled_index = index.iloc[:, [7, 6, 9, 10]].apply(zscore)
# 计算距离矩阵
dt = pdist(scaled_index, metric='euclidean')

#===================系统聚类===================
# 系统聚类
hc = linkage(dt, method='average')

# 展示聚类谱系图
plt.figure(figsize=(19, 11))
dendrogram(hc, labels=index.index, orientation='left', color_threshold=0)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Distance')
plt.ylabel('Amino Acids')
plt.show()

# 加载序列数据
S = list(SeqIO.parse("influenza_A_virus.fasta", "fasta"))

# 提取特征向量
def get_tz(sequences):
    lst = []
    for seq in sequences:
        n = np.zeros(6)
        weizhi = np.zeros(6)
        for i, s in enumerate(seq.seq):
            if s in "AILMVF":
                n[0] += 1
                weizhi[0] += i + 1
            if s == "Y":
                n[1] += 1
                weizhi[1] += i + 1
            if s in "QNSGTPW":
                n[2] += 1
                weizhi[2] += i + 1
            if s in "KRH":
                n[3] += 1
                weizhi[3] += i + 1
            if s in "ED":
                n[4] += 1
                weizhi[4] += i + 1
            if s == "C":
                n[5] += 1
                weizhi[5] += i + 1
        feature_vector = np.concatenate([n / len(seq), weizhi / n])
        lst.append(feature_vector)
    return np.array(lst)

# 运行函数所需时间
import time
start_time = time.time()
tz = get_tz(S)
end_time = time.time()
print(f"运行时间: {round(end_time - start_time,4)} 秒")

# 转换为数据框
tz_df = pd.DataFrame(tz, index=[rec.description for rec in S])

#===================相似性距离矩阵===================
# 计算距离矩阵之前进行标准化
scaled_tz = tz_df.apply(zscore)
# 计算距离矩阵
dt = pdist(scaled_tz, metric='euclidean')
dt_show = squareform(dt)

#===================系统聚类===================
# 系统聚类
hc = linkage(dt, method='average')

# 展示聚类谱系图
plt.figure(figsize=(19, 11))
dendrogram(
    hc,
    labels=tz_df.index,
    orientation='left',
    color_threshold=0,
    leaf_font_size=8  # 调整标签字体大小
)
# 设置边距
plt.subplots_adjust(left=0.1, right=0.2, top=0.2, bottom=0.1)  # 调整边距
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Linkage Distance')
plt.ylabel('Sequences')
plt.tight_layout()
plt.show()
