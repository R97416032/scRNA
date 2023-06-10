import csv
from utils.ungz import ungz
#解压文件
file_name="v43/basic/gencode.v43.basic.annotation.gtf.gz"

#读取gtf文件
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
                res.append(temp)
    #去除重复
    genecode=[]
    names=[]
    with open(path.replace(".gtf",".csv"), "w", newline='') as file:
        writer=csv.writer(file)
        cols=['type','gene_id','gene_type','gene_name']
        writer.writerow(cols)
        for r in res:
            if(r[-1] not in names):
                genecode.append(r)
                names.append(r[-1])
                writer.writerow(r)
ungz(file_name)



