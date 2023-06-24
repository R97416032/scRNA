import pandas as pd


def toWavalen(path):
    data=pd.read_csv(path)
    sample=data.columns[1:]
    geneid=data["gene_id"].values
    df=pd.DataFrame(columns=["gene_id",*sample])
    df["gene_id"]=geneid
    for s in sample:
        min_expression=min(data[s].values)
        max_expression=max(data[s].values)
        res=[]
        for v in data[s].values:
            res.append("%.2f"%((v-min_expression)*(780-380)/(max_expression-min_expression)+380))
        df[s]=res
    df.to_csv("wavelength.csv",index=None)
toWavalen("scFPKM.csv")