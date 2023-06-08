import csv
import numpy as np
import gzip
#解压文件
file_name="v43/basic/gencode.v43.basic.annotation.gtf.gz"
def ungz(file_name):
    f_name = file_name.replace(".gz", "")
    g_file = gzip.GzipFile(file_name)
    # 创建gzip对象
    open(f_name, "wb+").write(g_file.read())
    # gzip对象用read()打开后，写入open()建立的文件里。
    g_file.close()  # 关闭gzip对象
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



