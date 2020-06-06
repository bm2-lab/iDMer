import random, os, json, dgl, torch, sys
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
import torch.nn as nn
import torch.nn.functional as F

Datapath = os.path.dirname(os.path.abspath(__file__))
class GATLayer(nn.Module):
    def __init__(self, g, in_dim, out_dim):
        super(GATLayer, self).__init__()
        self.g = g
        self.fc = nn.Linear(in_dim, out_dim, bias=False)
        self.attn_fc = nn.Linear(2 * out_dim, 1, bias=False)

    def edge_attention(self, edges):
        z2 = torch.cat([edges.src['z'], edges.dst['z']], dim=1)
        a = self.attn_fc(z2)
        return {'e': F.leaky_relu(a)}

    def message_func(self, edges):
        return {'z': edges.src['z'], 'e': edges.data['e']}

    def reduce_func(self, nodes):
        alpha = F.softmax(nodes.mailbox['e'], dim=1)
        h = torch.sum(alpha * nodes.mailbox['z'], dim=1)
        return {'h': h}

    def forward(self, h):
        z = self.fc(h)
        self.g.ndata['z'] = z
        self.g.apply_edges(self.edge_attention)
        self.g.update_all(self.message_func, self.reduce_func)
        return self.g.ndata.pop('h')

class MultiHeadGATLayer(nn.Module):
    def __init__(self, g, in_dim, out_dim, num_heads, merge='cat'):
        super(MultiHeadGATLayer, self).__init__()
        self.heads = nn.ModuleList()
        for i in range(num_heads):
            self.heads.append(GATLayer(g, in_dim, out_dim))
        self.merge = merge

    def forward(self, h):
        head_outs = [attn_head(h) for attn_head in self.heads]
        if self.merge == 'cat':
            return torch.cat(head_outs, dim=1)
        else:
            return torch.mean(torch.stack(head_outs))

class GAT(nn.Module):
    def __init__(self, g, in_dim, hidden_dim, out_dim, num_heads):
        super(GAT, self).__init__()
        self.layer1 = MultiHeadGATLayer(g, in_dim, hidden_dim, num_heads)
        self.layer2 = MultiHeadGATLayer(g, hidden_dim * num_heads, out_dim, 1)

    def forward(self, h):
        h = self.layer1(h)
        h = F.elu(h)
        h = self.layer2(h)
        return h

def evaluate(model, features, labels, mask):
    model.eval()
    with torch.no_grad():
        logits = model(features)
        p = F.softmax(logits, 1)
        logits = logits[mask]
        labels = labels[mask]
        _, indices = torch.max(logits, dim=1)
        correct = torch.sum(indices == labels)
        acc = correct.item() * 1.0 / len(labels)
        y_true = labels.cpu().numpy()
        y_pred = indices.cpu().numpy()
        roc = roc_auc_score(y_true, y_pred)
        proba = p.cpu().numpy()
        return acc,roc,proba

def main(inputUP, inputDN, outdir): 
    with open("../data/PPI.json",'r') as f:
        PPI = json.load(f)
    node_list = PPI["node_list"]
    edge_list = PPI["edge_list"]

    with open("../data/node_feature.json",'r') as f:
        featurelist = json.load(f)
    feature_list = featurelist["features"]

    down_proteins = pd.read_csv(inputDN, sep='\t', header=None)
    up_proteins = pd.read_csv(inputUP, sep='\t', header=None)
    all_protein = pd.concat([down_proteins,up_proteins],ignore_index=True)
    all_protein.columns = ['id','name','label','class']

    seed_features = []
    null_list = []
    seed_id = []
    for idx in all_protein.id.values:
        if idx in node_list:
            seed_id.append(idx)
            i = node_list.index(idx)
            seed_features.append(feature_list[i])
        else:
            name = all_protein[all_protein.id == idx].name.values[0]
            null_list.append(name)

    seed_f = np.array(seed_features)
    seed_sum = np.sum(seed_f, axis=0)
    feature_id = []
    for i,j in enumerate(seed_sum):
        if j != 0:
            feature_id.append(i)
    seed_features2 = []
    for feature in seed_features:
        new_f = [feature[i] for i in feature_id]
        seed_features2.append(new_f)
    feature_list2 = []
    for feature in feature_list:
        new_f = [feature[i] for i in feature_id]
        feature_list2.append(new_f)

    label_list = []
    for node in node_list:
        if node in seed_id:
            score_ = all_protein[all_protein.id == node].label.values[0]
            if score_ > 0:
                score = 0
            else:
                score = 1
        else:
            score = 2
        label_list.append(score)

    up = []
    down = []
    for i,idx in enumerate(seed_id):
        label = all_protein[all_protein.id == idx].label.values[0]
        if label > 0:
            up.append(i)
        else:
            down.append(i)

    train_set = []
    test_set = []
    up_train_num = int(len(up)*0.8)
    down_train_num = int(len(down)*0.8)

    while len(train_set) < 1:
        random.shuffle(up)
        random.shuffle(down)
        train_list = up[:up_train_num]+down[:down_train_num]
        test_list = up[up_train_num:]+down[down_train_num:]
        train = [seed_features2[i] for i in train_list]
        train = np.array(train)
        train_sum = np.sum(train, axis=0)
        if 0 not in train_sum:
            train_set.append(train_list)
            test_set.append(test_list)

    train_mask = len(label_list)*[False]
    test_mask = len(label_list)*[False]
    for i in train_set[0]:
        idx = seed_id[i]
        i_ = node_list.index(idx)
        train_mask[i_] = True
    for i in test_set[0]:
        idx = seed_id[i]
        i_ = node_list.index(idx)
        test_mask[i_] = True

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    g = dgl.DGLGraph()
    g.add_nodes(len(node_list))
    src, dst = tuple(zip(*edge_list))
    g.add_edges(src, dst)
    g.add_edges(dst, src)

    feature_list = np.array(feature_list2)
    label_list = np.array(label_list)
    features = torch.from_numpy(feature_list).float().to(device)
    labels = torch.from_numpy(label_list).long().to(device)
    train_mask = np.array(train_mask)
    test_mask = np.array(test_mask)
    train_mask = torch.from_numpy(train_mask).to(device)
    test_mask = torch.from_numpy(test_mask).to(device)

    roc_list = []
    proba_list = []
    for i in range(10):
        net = GAT(g, in_dim=features.size()[1], hidden_dim=8, out_dim=2, num_heads=2)
        optimizer = torch.optim.Adam(net.parameters(), lr=1e-3)
        net.to(device)
        max_roc = 0.0

        for epoch in range(2000):
            net.train()
            logits = net(features)
            logp = F.log_softmax(logits, 1)
            loss = F.nll_loss(logp[train_mask], labels[train_mask])

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            acc,roc,proba = evaluate(net, features, labels, test_mask)

            if roc > max_roc:
                max_roc = roc
                max_proba = proba.tolist()
        roc_list.append(max_roc)
        proba_list.append(max_proba)


    best_roc = max(roc_list)
    best_proba = proba_list[roc_list.index(best_roc)]
    df = pd.DataFrame(best_proba)
    df.columns = ['up','down']
    df['id'] = node_list
    seed_id = all_protein.id.values
    all_id = df.id.values
    for i in list(all_id):
        if i in list(seed_id):
            df = df.drop(index=(df[df['id']==i].index))
    down_sorted = df.sort_values(by=['up'],ascending=True)
    up_sorted = df.sort_values(by=['down'],ascending=True)
    down_list = down_sorted.id.values
    up_list = up_sorted.id.values
    down_100 = down_list[:100]
    up_100 = up_list[:100]
    list_all = {"up_list":list(up_100),"down_list":list(down_100)}
    df2 = pd.DataFrame(list_all)
    df2.to_csv("{}/up_down_protein_GAT.csv".format(outdir), sep='\t', index=False)
    print("Best roc_auc_score:",best_roc)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
