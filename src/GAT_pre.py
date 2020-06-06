import numpy as np
import pandas as pd
import json, os

Datapath = os.path.dirname(os.path.abspath(__file__))
def main():
    
    links = []
    filein = '{}/data/9606.protein.links.v11.0.txt'.format(Datapath)
    with open(filein, "r") as f:
        f.readline()
        data = f.readlines()
        for data_ in data:
            link = data_.split()
            score = int(link[2])
            if score >= 900:
                links.append(link)

    id_name = {}
    filein = '{}/data/9606.protein.info.v11.0.txt'.format(Datapath)
    with open(filein, "r") as f:
        f.readline()
        data = f.readlines()
        for data_ in data:
            protein = data_.split()
            id_name.update({protein[0]:protein[1]})

    nodes = []
    for link in links:
        nodes += link[:2]
    node_list = list(set(nodes))

    edge_list = []
    for link in links:
        edge = (node_list.index(link[0]),node_list.index(link[1]))
        edge_list.append(edge)

    name_list = [id_name[node] for node in node_list]

    filein = '{}/data/gene2entrez.tsv'.format(Datapath)
    gene2entrez = pd.read_csv(filein, sep='\t')

    def get_type_list(GO_dict):
        type_list = []
        for name in name_list:
            if len(gene2entrez[gene2entrez.SYMBOL == name]) == 0:
                type_ = []
            else:
                ID = gene2entrez[gene2entrez.SYMBOL == name].ENTREZID.values[0]
                type_ = []
                for k in GO_dict:
                    if str(ID) in GO_dict[k]:
                        type_.append(k)
            type_list.append(type_)
        return type_list

    def get_GO(filename):
        dict_ = {}
        name = []
        with open(filename, "r") as f:
            f.readline()
            data = f.readlines()
            for data_ in data:
                l = data_.split('\t')
                if l[2] == '0':
                    continue
                l_ = l[-1].split("/")
                dict_.update({l[0]:l_})
                name.append(l[0])
            list_ = get_type_list(dict_)
        return list_,name
    
    MF_list, MF_name = get_GO('{}/data/MF_level2.txt'.format(Datapath))
    BP_list, BP_name = get_GO('{}/data/CC_level2.txt'.format(Datapath))
    CC_list, CC_name = get_GO('{}/data/BP_level2.txt'.format(Datapath))
    GO_name = MF_name+BP_name+CC_name
    feature_list = []
    for i in range(len(node_list)):
        feature = len(GO_name)*[0]
        GO_MF = [GO_name.index(GO) for GO in MF_list[i]]
        GO_BP = [GO_name.index(GO) for GO in BP_list[i]]
        GO_CC = [GO_name.index(GO) for GO in CC_list[i]]
        GO = GO_MF+GO_BP+GO_CC
        for j in GO:
            feature[j] = 1
        feature_list.append(feature)

    PPI = {"node_list":node_list,"edge_list":edge_list}
    PPI = json.dumps(PPI)
    fileout = '{}/data/PPI.json'.format(Datapath)
    with open(fileout,'w') as f:
        f.write(PPI)

    node_feature = {"features":feature_list}
    node_feature = json.dumps(node_feature)
    fileout = '{}/data/node_feature.json'.format(Datapath)
    with open(fileout,'w') as f:
        f.write(node_feature)

if __name__ == "__main__":
    main()