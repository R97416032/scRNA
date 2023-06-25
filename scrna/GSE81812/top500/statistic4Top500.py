import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from utils.showvalue import show_values


def getTop500(path,key="down"):
    data=pd.read_csv(path)
    if(key=="up"):
        data=data.sort_values(by="counts",ascending=False)
    data=data.head(500)
    gene_type=set(data["gene_type"].tolist())
    gene_type=[*gene_type]
    gene_type.sort()
    res=[]
    avg=[]
    for gt in gene_type:
        res.append(sum(data[data["gene_type"] == gt]["counts"].values))
        avg.append(sum(data[data["gene_type"]==gt]["counts"].values)/500)
    df=pd.DataFrame()
    df["gene_type"]=gene_type
    df["res"]=res
    df["avg"]=avg
    df=df.sort_values(by="res",ascending=False)
    return df["gene_type"].values,df["res"].values,df["avg"].values




path="../GSE81812_Normalized_counts_eachgene_result.csv"
#绘制图像最少的前500
gene_type,res,avg=getTop500(path)
ax0=sns.barplot(x=gene_type,y=res)
ax0.set_xticklabels(labels=gene_type,rotation = 90,fontsize = 8)
show_values(ax0)
plt.title("后500counts的各个类型基因分布")
plt.tight_layout()
plt.show()
# plt.savefig("../../../pic/00最后500counts的各个类型基因分布.jpg")
plt.close()
# ax0=sns.barplot(x=data["gene_id"].values,y=data["counts"].values)
# ax0.set_xticklabels(labels=data["gene_id"].values,rotation = 90,fontsize = 8)
# show_values(ax0)
# plt.title("后500各个基因的counts")
# plt.tight_layout()
# plt.show()
# plt.close()

##绘制饼图
ex=[]
colors=sns.color_palette("pastel")
i=0.1
for r in res:
    if r/sum(res) <0.1:
        ex.append(i)
        i+=0.04
    else:
        ex.append(0)
patches,l_text,p_text=plt.pie(res,labels=gene_type,colors=colors,autopct = '%1.2f%%',
        pctdistance = 1.1,labeldistance = 1.2,explode=ex)
for t in p_text:
    t.set_size(7)

for t in l_text:
    t.set_size(8)
plt.tight_layout()
plt.title("后500counts的各个类型基因分布饼图")

plt.show()
# plt.savefig("../../../pic/00最后500counts的各个类型基因分布饼图.jpg",dpi=300)
plt.close()
#绘制图像最多的前500
gene_type,res,avg=getTop500(path,"up")
ax1=sns.barplot(x=gene_type,y=res)
ax1.set_xticklabels(labels=gene_type,rotation = 90,fontsize = 8)
show_values(ax1)
plt.tight_layout()
plt.title("前500counts的各个类型基因分布")
plt.show()
# plt.savefig("../../../pic/01前500counts的各个类型基因分布.jpg")
plt.close()

##饼图
ex=[]
colors=sns.color_palette("pastel")
i=0.1
for r in res:
    if r/sum(res) <0.1:
        ex.append(i)
        i+=0.04
    else:
        ex.append(0)
patches,l_text,p_text=plt.pie(res,labels=gene_type,colors=colors,autopct = '%1.2f%%',
        pctdistance = 1.1,labeldistance = 1.2,explode=ex)
for t in p_text:
    t.set_size(7)

for t in l_text:
    t.set_size(8)
plt.tight_layout()
plt.title("前500counts的各个类型基因分布饼图")
plt.show()


# plt.savefig("../../../pic/01前500counts的各个类型基因分布饼图.jpg")
plt.close()

