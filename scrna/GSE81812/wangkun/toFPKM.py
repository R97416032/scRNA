import csv
from utils.ungz import ungz
import numpy as np
import pandas as pd

def read(path):
    with open(path, 'r') as gtf:
        res=[]
        for g in gtf:
            if g[0]=='c':
                glist=g.strip().split('\t')
                temp = [glist[2]]
                for infos in glist[-1].split(';'):
                    info=infos.split()
                    if(len(info)==0):
                        continue
                    if(info[0]=='gene_id'):
                        temp.append(info[1].strip('"').split('.')[0])
                    elif(info[0]=='gene_type'):
                        temp.append(info[1].strip('"'))
                    elif(info[0]=='gene_name'):
                        temp.append(info[1].strip('"'))
                        continue
                temp.append(glist[3])
                temp.append(glist[4])
                res.append(temp)
    with open("gencode.v43.basic.annotation_withlen.csv", "w", newline='') as file:
        writer=csv.writer(file)
        cols=['type','gene_id','gene_type','gene_name','gene_start','gene_end']
        writer.writerow(cols)
        for r in res:
            writer.writerow(r)
def drop_duplicate():
    read("../../../genecode/v43/basic/gencode.v43.basic.annotation.gtf")

    data=pd.read_csv("gencode.v43.basic.annotation_withlen.csv")

    type=set(data["gene_type"].values)
    for t in type:
        if "pseudogene" in t:
                data.drop(index=data[data["gene_type"]==t].index,inplace=True)
    data.to_csv("gencode.v43.basic.annotation_withlen.csv",index=None)

    data=pd.read_csv("gencode.v43.basic.annotation_withlen.csv")
    ns=data["gene_name"].values.tolist()
    names=list(set(ns))
    names.sort(key=ns.index)

    for n in ns :
        if n not in names:
            names.append(n)
    print(names)
    for n in names:
        data.loc[data[data["gene_name"] == n].index, "gene_start"]=min(data.loc[data[data["gene_name"]==n].index,"gene_start"].values)
        data.loc[data[data["gene_name"] == n].index, "gene_end"]=max(data.loc[data[data["gene_name"]==n].index,"gene_end"].values)
        print(n)
    data.to_csv("gencode.v43.basic.annotation_withlen.csv",index=None)
    # 去除重复
    data.drop_duplicates(["gene_name"],inplace=True)
    data.to_csv("gencode.v43.basic.annotation_withlenv2.csv",index=None)

def drop_pseudo_and_un(path):
    data=pd.read_csv(path)
    type = set(data["gene_type"].values)
    for t in type:
        if "pseudogene" in t:
            data.drop(index=data[data["gene_type"] == t].index, inplace=True)
        elif  "unknow" in t:
            data.drop(index=data[data["gene_type"] == t].index, inplace=True)
    data.to_csv("scdata.csv", index=None)
# drop_pseudo_and_un("../GSE81812_Normalized_counts.csv")
def toFPKM(path):
    scdata=pd.read_csv(path)
    genecode=pd.read_csv("gencode.v43.basic.annotation_withlenv2.csv")
    genecode_id=genecode["gene_id"].values
    sample=scdata.columns.values[2:]
    geneid=scdata["gene_id"].values
    genelength=[]
    for id in geneid:
        if id in genecode_id:
            genelength.append(int(genecode.loc[genecode[genecode["gene_id"]==id].index,"gene_end"])-int(genecode.loc[genecode[genecode["gene_id"]==id].index,"gene_start"]))
        else:
            genelength.append(-1)
    data=pd.DataFrame(columns=["gene_id",*sample])

    data["gene_id"]=geneid
    for s in sample:
        fpkm=[]
        counts=scdata[s].values
        sumcounts=sum(counts)
        for c,l in zip(counts,genelength):#zip将两个列表合并
            if l !=-1:
                fpkm.append('%.2f'%(c*10e9/(l*sumcounts)))
            else:
                fpkm.append(-1)
        data[s]=fpkm
    data.to_csv("scFPKM.csv",index=None)
toFPKM("scdata.csv")
