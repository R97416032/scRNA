import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def fpkmtocolor(path):
    data = pd.read_csv(path)
    sample = data.columns[1:]
    geneid = data["gene_id"].values
    df = pd.DataFrame(columns=["gene_id", *sample])
    df["gene_id"] = geneid
    for s in sample:
        res = []
        for v in data[s].values:
            if v>=380 and v<420:
                res.append(7)
            elif v>=420 and v<450:
                res.append(6)
            elif v>=450 and v<490:
                res.append(5)
            elif v>=490 and v<560:
                res.append(4)
            elif v>=560 and v<590:
                res.append(3)
            elif v>=590 and v<620:
                res.append(2)
            elif v>=620 and v<=780:
                res.append(1)
            else:
                print(-1)
                res.append(-1)
        df[s] = res
    df.to_csv("color.csv", index=None)

def heatmap(path):
    data = pd.read_csv(path)
    sample = data.columns[1:]
    geneid = data["gene_id"].values
    df=pd.DataFrame(columns=["sample","gene_id","color"])
    samples=[]
    colors=[]
    geneids=[]
    for s in sample:
        for id,color in zip(geneid,data[s].values):
            samples.append(s)
            geneids.append(id)
            colors.append(color)
    df["sample"]=samples
    df["gene_id"]=geneids
    df["color"]=colors
    df.to_csv("colors.csv",index=None)
    df=pd.read_csv("colors.csv")
    sns.heatmap(df.pivot("sample","gene_id","color"),xticklabels=1,yticklabels=1)
    plt.show()
fpkmtocolor("wavelength.csv")
heatmap("color.csv")
