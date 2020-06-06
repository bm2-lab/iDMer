import argparse
parser = argparse.ArgumentParser(description='It is used for Compound Identification',formatter_class=argparse.RawDescriptionHelpFormatter,add_help=True)
subparsers = parser.add_subparsers(help='commands')

### experiment host factors
exp = subparsers.add_parser('exp',help='host factors collected from experiment',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
Req = exp.add_argument_group('Required')
Req.add_argument('-VDN', action='store', metavar='<VTP reliance gene>', required=True)
Req.add_argument('-VUP', action='store', metavar='<VTP restriction gene>', required=True)
Req.add_argument('-EDN', action='store', metavar='<EHF reliance gene>', required=True)
Req.add_argument('-EUP', action='store', metavar='<EHF restriction gene>', required=True)
Opt = exp.add_argument_group('Optional')
Opt.add_argument('-output', metavar='<output_dir>',action='store', default='iDMer', help='default iDMer')
Opt.add_argument('-name', metavar='<CMap_name>',action='store', default='CMap', help='default CMap')
Opt.add_argument('-g','--GAT',action='store_true',default=False, help='implement GAT, will take one day')
Opt.add_argument('-h','--help',action='help', help='show this help message and exit')

### in silico prediction host factors
denovo = subparsers.add_parser('denovo',help='host factors obtained from in silico',add_help=False,formatter_class=argparse.RawTextHelpFormatter)
Req = denovo.add_argument_group('Required')
Req.add_argument('-virus', action='store', metavar='<viral gene information in fasta format>', required=True)
Req.add_argument('-host', action='store', metavar='<viral target proteins in fasta format>', required=True)
Req.add_argument('-config', action='store', metavar='<file indicating the class of the predicted host factors>', required=True)
Req.add_argument('-EDN', action='store',metavar='<EHF reliance gene>',required=True)
Req.add_argument('-EUP', action='store',metavar='<EHF restriction gene>',required=True)
Opt = denovo.add_argument_group('Optional')
Opt.add_argument('-output',metavar='output_dir',action='store', default='iDMer', help='default iDMer')
Opt.add_argument('-name', metavar='<CMap_name>',action='store', default='CMap', help='default CMap')
Opt.add_argument('-VDN', metavar= '<VTP reliance gene>',action='store', default= 'VTPs_DN.tsv', help = 'default VTPs_DN.tsv')
Opt.add_argument('-VUP', metavar= '<VTP restriction gene>',action='store', default= 'VTPs_UP.tsv', help = 'default VTPs_UP.tsv')
Opt.add_argument('-g','--GAT',action='store_true',default=False, help='implement GAT, will take one day')
Opt.add_argument('-h','--help',action='help', help='show this help message and exit')

args = parser.parse_args()
import os, subprocess, time
import pandas as pd, numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

Datapath = os.path.dirname(os.path.abspath(__file__))
ensembl2geneName = {}
geneName2ensembl = {}

file = '{}/data/9606.protein.info.v11.0.txt'.format(Datapath)
with open(file, 'r') as fin:
    fin.readline()
    for line in fin:
        lines = line.strip().split('\t')
        ensembl2geneName[lines[0]] = lines[1]
        geneName2ensembl[lines[1]] = lines[0]

def HVPPI(args):
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    host_dict = {}
    virus_dict = {}
    with open(args.virus, 'r') as fin:
        for title, seq in SimpleFastaParser(fin):
            virus_dict[title] = seq
    with open(args.host, 'r') as fin:
        for title, seq in SimpleFastaParser(fin):
            host_dict[title] = seq
    fileout1 = '{}/HVPPI.fa'.format(args.output)
    fileout2 = '{}/HVPPI.tsv'.format(args.output)
    with open(fileout1,'w') as fout1, open(fileout2, 'w') as fout2:
        for i in host_dict:
            fout1.write('{}\t{}\n'.format(i, host_dict[i]))
        for i in virus_dict:
            fout1.write('{}\t{}\n'.format(i, virus_dict[i]))
        for i in host_dict:
            for j in virus_dict:
                fout2.write('{}\t{}\n'.format(i, j))
    cmd = 'python /home/src/doc2vec_rf.py  {}  {}  0.9  {}'.format(fileout2, fileout1, args.output)
    subprocess.call(cmd, shell=True)
    filein = '{}/PPI_prediction_result.out'.format(args.output)
    prediction = pd.read_csv(filein, sep='\t')
    dat = pd.DataFrame(prediction.groupby('Pro1ID')['Score'].max())
    dat['Pro1ID'] = dat.index
    dat.sort_values(by='Score', ascending=False, inplace=True)
    prediction_genelist = dat['Pro1ID'][:150].tolist()
    config_dict = {}
    with open(args.config, 'r') as fin:
        for line in fin:
            lines = line.strip().split('\t')
            config_dict[lines[0]] = int(lines[-1])
    fileDN = '{}'.format(args.VDN)
    fileUP = '{}'.format(args.VUP)
    with open(fileDN, 'w') as fout1, open(fileUP, 'w') as fout2:
        for i in prediction_genelist:
            if i in config_dict and config_dict[i] == -1:
                fout1.write('{}\n'.format(i))
            if i in config_dict and config_dict[i] == 1:
                fout2.write('{}\n'.format(i))

def getlist(args):
    total_list = []
    down_proteins = []
    up_proteins = []
    with open(args.VDN, 'r') as fin:
        for line in fin:
            protein = line.strip()
            if protein in geneName2ensembl and protein not in total_list:
                ensembl = geneName2ensembl[protein]
                down_proteins.append([ensembl, protein, -1, 'VTPs'])
                total_list.append(protein)
    with open(args.EDN, 'r') as fin:
        for line in fin:
            protein = line.strip()
            if protein in geneName2ensembl and protein not in total_list:
                ensembl = geneName2ensembl[protein]
                down_proteins.append([ensembl, protein, -1, 'EHFs'])
                total_list.append(protein)

    with open(args.VUP, 'r') as fin:
        for line in fin:
            protein = line.strip()
            if protein in geneName2ensembl and protein not in total_list:
                ensembl = geneName2ensembl[protein]
                up_proteins.append([ensembl, protein, 1, 'VTPs'])
                total_list.append(protein)
    with open(args.EUP, 'r') as fin:
        for line in fin:
            protein = line.strip()
            if protein in geneName2ensembl and protein not in total_list:
                ensembl = geneName2ensembl[protein]
                up_proteins.append([ensembl, protein, 1, 'EHFs'])
                total_list.append(protein)
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    fileout = '{}/down_proteins.tsv'.format(args.output)
    with open(fileout, 'w') as fout:
        for i in down_proteins:
            fout.write('{}\t{}\t{}\t{}\n'.format(i[0], i[1], i[2], i[3]))

    fileout = '{}/up_proteins.tsv'.format(args.output)
    with open(fileout, 'w') as fout:
        for i in up_proteins:
            fout.write('{}\t{}\t{}\t{}\n'.format(i[0], i[1], i[2], i[3]))

def getGATlist():
    count1 = len(open('{}/up_proteins.tsv'.format(args.output),'r').readlines())
    count2 = len(open('{}/down_proteins.tsv'.format(args.output),'r').readlines())
    if count1 + count2 >= 100:
        GAT = 9
    else:
        GAT = 4
    filein = '{}/up_down_protein_GAT.csv'.format(args.output)
    fileup = '{}/up_proteins.GAT.tsv'.format(args.output)
    filedn = '{}/down_proteins.GAT.tsv'.format(args.output)
    with open(filein, 'r') as fin, open(fileup, 'w') as fup, open(filedn, 'w') as fdn:
        fin.readline()
        up_GAT = 0
        down_GAT = 0
        for line in fin:
            lines = line.strip().split('\t')
            ensembl = lines[0]
            if ensembl in ensembl2geneName:
                protein = ensembl2geneName[ensembl]
                if up_GAT <= GAT:
                    fup.write('{}\t{}\t{}\t{}\n'.format(ensembl, protein, 1, 'GAT'))
                    up_GAT += 1
            ensembl = lines[1]   ### down
            if ensembl in ensembl2geneName:
                protein = ensembl2geneName[ensembl]
                if down_GAT <= GAT:
                    fdn.write('{}\t{}\t{}\t{}\n'.format(ensembl, protein, -1, 'GAT'))
                    down_GAT += 1
    cmd = 'cat {} >>  {}/up_proteins.tsv'.format(fileup, args.output)
    subprocess.call(cmd, shell=True)
    cmd = 'cat {} >>  {}/down_proteins.tsv'.format(filedn, args.output)
    subprocess.call(cmd, shell=True)

def doGAT(args):
    print ('\n********usually take one day to complete*********\n')
    PPI = '{}/data/PPI.json'.format(Datapath)
    node = '{}/data/node_feature.json'.format(Datapath)
    if not os.path.isfile(PPI) or not os.path.isfile(node):
        cmd = 'python {}/src/GAT_pre.py'.format(Datapath)
        print (cmd)
        subprocess.call(cmd, shell=True)
    cmd = 'python {}/src/GAT_train.py {}/up_proteins.tsv  {}/down_proteins.tsv {}'.format(Datapath ,args.output, args.output, args.output)
    print (cmd)
    subprocess.call(cmd, shell=True)
    getGATlist()

def postCMap(name):
    print ('\n********usually take half an hour to complete********\n')
    gene2entrez = {}
    filein = '{}/data/Broad_LINCS_gene_info.txt'.format(Datapath)
    with open(filein, 'r') as fin:
        fin.readline()
        for line in fin:
            lines = line.strip().split('\t')
            gene2entrez[lines[1]] = lines[0]
    temp = pd.read_table(filein)
    bing = temp['pr_gene_symbol'][temp['pr_is_bing'] == 1].tolist()
    up_filein = '{}/up_proteins.tsv'.format(args.output)
    up_fileout = '{}/uptag.gmt'.format(args.output)
    with open(up_filein, 'r') as fin, open(up_fileout, 'w') as fout:
        up_tag = [line.strip().split('\t')[1] for line in fin]
        up_tag = [i for i in up_tag if i in bing]
        up_tag = [gene2entrez[i] for i in up_tag if i in gene2entrez]
        if len(up_tag) >= 150:
            up_tag = up_tag[:150]
        fout.write('TAG_UP\t\t{}\n'.format('\t'.join(up_tag)))
        print('up_tag has {} protein\n'.format(len(up_tag)))

    dn_filein = '{}/down_proteins.tsv'.format(args.output)
    dn_fileout = '{}/dntag.gmt'.format(args.output)
    with open(dn_filein, 'r') as fin, open(dn_fileout, 'w') as fout:
        dn_tag = [line.strip().split('\t')[1] for line in fin]
        dn_tag = [i for i in dn_tag if i in bing]
        dn_tag = [gene2entrez[i] for i in dn_tag if i in gene2entrez]
        if len(dn_tag) >= 150:
            dn_tag = dn_tag[:150]
        fout.write('TAG_DN\t\t{}\n'.format('\t'.join(dn_tag)))
        print('dn_tag has {} protein\n'.format(len(dn_tag)))
    cmd = 'curl -i -X POST -H "user_key: 26dc7a7f528b22f0023dbef6311db6ef" -H "Content-Type: multipart/form-data"  ' \
          '-F "tool_id=sig_gutc_tool"  -F "uptag-cmapfile=@./{}"  -F "name={}"   ' \
          '-F "dntag-cmapfile=@./{}"  -F "data_type=L1000"  -F "dataset=Touchstone" ' \
          '"https://api.clue.io/api/jobs" '.format(up_fileout, name, dn_fileout)
    result = subprocess.getoutput(cmd)
    lines = result.strip().split()[-1].split(',')[-2].split('/')
    time = lines[-3]
    id = lines[-2]
    link = 'https://s3.amazonaws.com/data.clue.io/api/1632738@tongji.edu.cn/results/{}/{}/{}.tar.gz'.format(time, id, id)
    cmd = 'wget -c {} -O {}/{}.tar.gz'.format(link, args.output, name)
    with open('{}/CMap_link.txt'.format(args.output), 'w') as fout:
        fout.write('{}\n'.format(cmd))


def for_getdrug(name):
    filein = '{}/{}/matrices/gutc/ps_pert_summary.gctx'.format(args.output, name)
    fileout = '{}/ps_pert_summary.gct'.format(args.output)
    cmd = '/root/miniconda2/bin/gctx2gct -filename {} -output_filepath  ./{}'.format(filein, fileout)
    subprocess.call(cmd, shell=True)
    dat = pd.read_table(fileout, skiprows=2)
    dat.index = dat.pert_iname
    dat = dat[dat['pert_type'] == 'trt_cp']
    dat = dat[~dat['pert_iname'].str.startswith('BRD-')]
    dat.sort_values(by='TAG', inplace=True, ascending=False)
    dat = dat.iloc[:,[1,2,4]]
    dat.to_csv('{}/CMap_raw_result.tsv'.format(args.output), sep='\t', header=True, index=False)

def getdrug(args):
    with open('{}/CMap_link.txt'.format(args.output), 'r') as fin:
        cmd = fin.readline().strip()
        subprocess.call(cmd, shell=True)
    id = cmd.split()[2].split('/')[-2]
    name = args.name
    if not os.path.isdir(name):
        cmd = 'tar -zxvf  {}/{}.tar.gz'.format(args.output, name)
        subprocess.call(cmd, shell=True)
        cmd = 'mv {} {}/{}'.format(id, args.output, name)
        subprocess.call(cmd, shell=True)
    for_getdrug(name)


def toxicity(args):
    filein1 = '{}/data/BRD2Toxicity.tsv'.format(Datapath)
    filein2 = '{}/CMap_raw_result.tsv'.format(args.output)
    fileout = '{}/CMap_tox.tsv'.format(args.output)
    dat1 = pd.read_csv(filein1, sep='\t')
    dat = pd.read_csv(filein2, sep='\t')
    dat.drop(labels=['pert_iname'], axis=1, inplace=True)
    tmp = pd.merge(left=dat, left_on='pert_id', right=dat1, right_on='pert_id')
    tmp.to_csv(fileout, sep='\t', index=False)

def isDenovo(args):
    try:
        if args.virus:
            return True
    except:
        return False

if __name__ == '__main__':
    if isDenovo(args):
        HVPPI(args)
    getlist(args)
    if args.GAT is True:
        doGAT(args)
    postCMap(args.name)
    time.sleep(1800)   ### usually take half an hour to complete compound identification
    getdrug(args)
    toxicity(args)
