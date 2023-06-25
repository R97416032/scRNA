import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from utils.showvalue import show_values


def get_entropy(path):
    data=pd.read_csv(path)
    gene_id=data["gene_id"]
    sample=data.columns[1:]
    res=[]
    for i in range(len(gene_id)):
        # print(sum(data.iloc[i,:].values[1:]==1))
        probability=[]
        entropy=0
        for j in range(0,8):
            if (sum(data.iloc[i,:].values[1:]==j))>0:
                entropy=entropy-(sum(data.iloc[i,:].values[1:]==j)/len(sample))*math.log2(sum(data.iloc[i,:].values[1:]==j)/len(sample))
        res.append(entropy)
    df=pd.DataFrame(columns=["gene_id","entropy"])
    df["gene_id"]=gene_id
    df["entropy"]=res
    df.sort_values(by="entropy",inplace=True,ascending=False)
    df.to_csv("entropy.csv",index=None)

def get_type(path):
    data=pd.read_csv(path)
    genecode=pd.read_csv("gencode.v43.basic.annotation_withlenv2.csv")
    data=data[data["entropy"].values>0]
    gene_id=data["gene_id"].values
    gene_type=[]
    for id in gene_id:
        gene_type.append(genecode[genecode["gene_id"]==id]["gene_type"].values[0])
    df=pd.DataFrame(columns=["gene_id","gene_type","entropy"])
    df["gene_id"]=gene_id
    df["gene_type"]=gene_type
    df["entropy"]=data["entropy"].values
    df.to_csv("entropy_withtype.csv",index=None)
def show(path):
    data=pd.read_csv(path)
    type=set(data["gene_type"].values)
    type=[*type]
    num=[]
    for t in type:
        num.append(sum(data["gene_type"].values==t))
    ax=sns.barplot(x=type,y=num)
    show_values(ax)
    plt.show()
# get_entropy("color.csv")
# get_type("entropy.csv")
show("entropy_withtype.csv")