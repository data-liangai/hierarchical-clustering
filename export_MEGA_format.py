import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from Bio import SeqIO
import time
import os

os.chdir("D:\hierarchical-clustering")
# 加载序列数据
sequences = list(SeqIO.parse("influenza_A_virus.fasta", "fasta"))

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
start_time = time.time()
tz = get_tz(sequences)
end_time = time.time()
print(f"运行时间: {round(end_time - start_time,4)} 秒")

# 转换为数据框
tz_df = pd.DataFrame(tz)

# 第一步：导出抬头格式
with open("dist.meg", "w") as f:
    f.write("#mega\n!Title:;\n")

# 第二步：导出序列描述
description = [f"#{rec.description.replace(' ', '_')}" for rec in sequences]
r_name = [f"[{i+1}]" for i in range(len(sequences))]

description_df = pd.DataFrame(description, index=r_name)
description_df.to_csv("dist.meg", mode="a", sep="\t", header=False, index=True)

# 第三步：导出矩阵的列抬头
c_name = '\t'.join(str(i+1) for i in range(len(sequences)))
with open("dist.meg", "a") as f:
    f.write(f"[{c_name}]\n")

tz_df.index = r_name

# 用欧氏距离定义样本间的距离
scaled_tz = (tz_df - tz_df.mean()) / tz_df.std()  # 标准化数据
dt = pdist(scaled_tz, metric='euclidean')
dt_matrix = squareform(dt)
# 将上三角矩阵和主对角线都赋值为NA
tri_upper_indices = np.triu_indices_from(dt_matrix)
dt_matrix[tri_upper_indices] = np.nan

# 第四步：导出矩阵
dt_df = pd.DataFrame(dt_matrix, index=r_name, columns=r_name)
dt_df.to_csv("dist.meg", mode="a", sep="\t", header=False, index=True, float_format="%.4f", na_rep="")
