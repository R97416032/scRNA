import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator


def reshapeData(path):
    data = pd.read_csv(path)
    sample = data.columns.values[1:]
    gene_type = data["gene_type"].values
    X = []
    Y = []
    res=[]
    for s in sample:
        for gt in gene_type:
            Y.append(gt)
            X.append(s)
        for c in data[s].values.tolist():
            res.append(float(c)/1000)
    df=pd.DataFrame()
    df["x"]=X
    df["Y"]=Y
    df["counts"]=res
    return X,res,Y

path="GSE81812_Normalized_counts_result.csv"
data = pd.read_csv(path)
sample = data.columns.values[1:]
gene_type = data["gene_type"].values
x,y,g=reshapeData(path)
print(y)
print(len(x),len(y),len(g))
plt.scatter(x, g, s=y, alpha=0.6)
plt.tick_params(axis='x', labelsize=4)
plt.xticks(rotation=90)
plt.show()
plt.close()

