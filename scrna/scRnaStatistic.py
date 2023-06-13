import pandas as pd


def scRnaStatistic(path):
    data=pd.read_csv(path)
    gene_type=set(data["gene_type"].values.tolist())
    gene_type=[*gene_type]
    gene_type.sort()
    sample=data.columns[1:]
    df=pd.DataFrame(columns=sample)
    df["gene_type"]=gene_type

    for sam in sample[1:]:
        res=[]
        for gt in gene_type:
            temp=data[data["gene_type"]==gt]
            res.append(sum(temp[sam].values))
        df[sam]=res
    df.to_csv(path.replace(".csv","_result.csv"),index=None)
    df=pd.DataFrame()
    gene_id=data["gene_id"].values.tolist()
    gene_type = data["gene_type"].values.tolist()
    df["gene_id"] = gene_id
    df["gene_type"] = gene_type
    res=[]
    for gid in gene_id:
        res.append(sum(data[data["gene_id"]==gid].values[0][2:]))
    df["counts"]=res
    res = []
    for gid in gene_id:
        res.append("%.4f%%"%(sum(data[data["gene_id"] == gid].values[0][2:])/sum(df["counts"])))
    df["ratio"] = res
    df = df.sort_values(by="counts")
    df.to_csv(path.replace(".csv","_eachgene_result.csv"),index=None)


def computePercentage(path):
    data=pd.read_csv(path)
    sample=data.columns
    df=pd.DataFrame(columns=sample)
    gene_type =data["gene_type"].values
    df["gene_type"]=gene_type
    for sam in sample[1:]:
        # print(sam)
        sumcount=sum(data[sam].values)
        res=[]
        for gt in gene_type:
            res.append("%.5f%%"%((data[data["gene_type"]==gt][sam].values/sumcount)[0]*100))
        df[sam]=res

    df.to_csv(path.replace("_counts_result.csv","_precentage_result.csv"))


path="GSE81812/GSE81812_Normalized_counts.csv"
result_path="GSE81812/GSE81812_Normalized_counts_result.csv"
scRnaStatistic(path)
computePercentage(result_path)