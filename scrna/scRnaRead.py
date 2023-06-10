import csv
from utils.ungz import ungz
import numpy as np
import pandas as pd
filename="GSE81812/GSE81812_Normalized_counts.txt.gz"
ungz(filename)
data=[]
with open("GSE81812/GSE81812_Normalized_counts.txt",'r') as file:
    for f in file:
        row=f.strip().replace('"',"").split(",")
        if(row[0]==''):
            row[0]="gene_id"
        else:
            row[0]=row[0].split(".")[0]
        data.append(row)
with open(filename.replace(".txt.gz",".csv"), "w", newline='') as file:
    writer=csv.writer(file)
    for d in data:
        writer.writerow(d)

data=pd.read_csv(filename.replace(".txt.gz",".csv"))
genecode=pd.read_csv("../genecode/v43/basic/gencode.v43.basic.annotation.csv")
gene_id=genecode["gene_id"].values.tolist()
gene_type=genecode["gene_type"].values.tolist()
gty=[]
for gid in data["gene_id"].values:
    if gid in gene_id:
        gty.append(gene_type[gene_id.index(gid)])
    else:
        gty.append("unknow")
data.insert(1,"gene_type",gty)
data.to_csv(filename.replace(".txt.gz",".csv"),index=None)