import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
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

def show_values(axs, orient="v", space=.01):
    def _single(ax):
        if orient == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height() + (p.get_height()*0.01)
                value = '{:.0f}'.format(p.get_height())
                ax.text(_x, _y, value, ha="center")
        elif orient == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height() - (p.get_height()*0.5)
                value = '{:.0f}'.format(p.get_width())
                ax.text(_x, _y, value, ha="left")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _single(ax)
    else:
        _single(axs)

fig = plt.figure()
#绘制图像最少的前500
gene_type,res,avg=getTop500(path)
ax0=sns.barplot(x=gene_type,y=res)
ax0.set_xticklabels(labels=gene_type,rotation = 90,fontsize = 8)
show_values(ax0)
plt.tight_layout()
plt.savefig("../../../pic/00最后500counts的各个类型基因分布.jpg",dpi=640)
plt.close()
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
    t.set_size(5)

for t in l_text:
    t.set_size(6)
plt.savefig("../../../pic/00最后500counts的各个类型基因分布饼图.jpg",dpi=640)
plt.close()
#绘制图像最多的前500
gene_type,res,avg=getTop500(path,"up")
ax1=sns.barplot(x=gene_type,y=res)
ax1.set_xticklabels(labels=gene_type,rotation = 90,fontsize = 8)
show_values(ax1)
plt.tight_layout()
plt.savefig("../../../pic/01前500counts的各个类型基因分布.jpg")
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
# for t in p_text:
#     t.set_size(5)
#
# for t in l_text:
#     t.set_size(6)
plt.savefig("../../../pic/01前500counts的各个类型基因分布饼图.jpg")
plt.close()