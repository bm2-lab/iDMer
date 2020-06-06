#-*-coding:utf-8-*-
import pandas as pd
import os, sys, warnings
warnings.filterwarnings('ignore')

Datapath = os.path.dirname(os.path.abspath(__file__))
BRD2SMILE = {} 
file = '{}/data/BRD2SMILE.tsv'.format(Datapath) 
with open(file, 'r') as fin:
    for line in fin:
        lines = line.strip().split('\t')
        BRD2SMILE[lines[0]] = lines[1]
      
def combine_rank(filein, output):
    if not os.path.isdir(output):
        os.makedirs(output)
    df_1 = pd.read_csv(filein, sep="\t", header= 0)
    df_2 = pd.read_csv("{}/data/CRS_tox.tsv".format(Datapath), sep="\t")
    df_1_10 = df_1[df_1["Toxicity"] == "NO_Toxicity"].iloc[0:10, ]
    df_2_10 = df_2[df_2["Toxicity"] == "NO_Toxicity"].iloc[0:10, ]
    pert_id_1 = list(df_1_10["pert_id"])
    pert_id_2 = list(df_2_10["pert_id"])
    pert_id_same=[]    #the same drugs
    pert_id_differ=[]  #the different drugs
    pert_id_smile=[]
    for id_1 in pert_id_1:
        for id_2 in pert_id_2:
            if id_1==id_2:
                pert_id_same.append([id_1,id_2])
            else:
                id_1_smile=BRD2SMILE[id_1]
                id_2_smile=BRD2SMILE[id_2]
                pert_id_smile.append([id_1,id_1_smile,id_2,id_2_smile])
                pert_id_differ.append([id_1,id_2])
    df_smile_info = pd.DataFrame(pert_id_smile)
    df_smile_info.to_csv("{}/compoundCombinationSmile.tsv".format(output), sep="\t", index=False, header=None)
    #system deepddi cmd
    cmd_1="python /home/deepddi/run_DeepDDI.py -i  {}/compoundCombinationSmile.tsv -o ./combination_tmp/".format(output)
    os.system(cmd_1)
    #analyze the result
    DDI_anta_types=["DDI type 18","DDI type 19","DDI type 20","DDI type 21","DDI type 22",
              "DDI type 23","DDI type 25","DDI type 26","DDI type 27","DDI type 28",
              "DDI type 29","DDI type 30","DDI type 31","DDI type 32"]
    df_drugs_types = pd.read_csv("./combination_tmp/Final_DDI_result.txt", sep="\t")
    list_last_result=[]
    for [id_1,id_2] in pert_id_differ:
        score_1=list(df_1_10[df_1_10["pert_id"]==id_1]["TAG"])[0]
        score_2=list(df_2_10[df_2_10["pert_id"]==id_2]["TAG"])[0]
        name_1=list(df_1_10[df_1_10["pert_id"]==id_1]["pert_iname"])[0]
        name_2 = list(df_2_10[df_2_10["pert_id"] == id_2]["pert_iname"])[0]
        pair_1=str(id_1)+"_"+str(id_2)
        pair_2=str(id_2)+"_"+str(id_1)
        df_type_1=df_drugs_types[df_drugs_types["Drug pair"]==pair_1]
        df_type_2=df_drugs_types[df_drugs_types["Drug pair"]==pair_2]
        if len(df_type_1)==0 and len(df_type_2)==0:
            type_1="none"
            type_2="none"
            antangonism="no"
            score=(score_1+score_2)/2
        if len(df_type_1)==1 and len(df_type_2)==1:
            type_1=list(df_type_1["DDI type"])[0]
            type_2=list(df_type_2["DDI type"])[0]
            sentence_1=list(df_type_1["Sentence"])[0]
            sentence_2=list(df_type_2["Sentence"])[0]
            if (type_1 not in DDI_anta_types) and (type_2 not in DDI_anta_types):
                antangonism = "no"
                score = (score_1 + score_2) / 2
            else:
                antangonism = "yes"
                score = 0
        if len(df_type_1)==1 and len(df_type_2)==0:
            type_1 = list(df_type_1["DDI type"])[0]
            sentence_1 = list(df_type_1["Sentence"])[0]
            type_2="none"
            sentence_2="none"
            if type_1 not in DDI_anta_types:
                antangonism = "no"
                score = (score_1 + score_2) / 2
            else:
                antangonism = "yes"
                score = 0
        if len(df_type_1) == 0 and len(df_type_2) == 1:
            type_2 = list(df_type_2["DDI type"])[0]
            sentence_2 = list(df_type_2["Sentence"])[0]
            type_1 = "none"
            sentence_1 = "none"
            if type_2 not in DDI_anta_types:
                antangonism = "no"
                score = (score_1 + score_2) / 2
            else:
                antangonism = "yes"
                score = 0
        list_last_result.append([id_1,name_1,id_2,name_2,type_1,sentence_1,type_2,sentence_2,antangonism,score])

    for [id_1,id_2] in pert_id_same:
        name_1 = list(df_1_10[df_1_10["pert_id"] == id_1]["pert_iname"])[0]
        name_2 = list(df_2_10[df_2_10["pert_id"] == id_2]["pert_iname"])[0]
        type_1 = "none"
        sentence_1 = "none"
        antangonism = "no"
        score=100
        list_last_result.append(
            [id_1, name_1, id_2, name_2, type_1, sentence_1, type_2, sentence_2, antangonism, score])

    result = pd.DataFrame(list_last_result, columns=["pert_id_1", "pert_iname_1","pert_id_2","pert_iname_2", "ddi_type_1", "sentence_1", "ddi_type_2", "sentence_2",
                                             "antangonism","combine_score"])
    result = result.sort_values("combine_score", ascending=False)
    fileout = '{}/compoundCombination.tsv'.format(output)
    result.to_csv(fileout, sep="\t", index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print ('\nError\npython {}  virus_results  output_dirName\n'.format(sys.argv[0]))
        sys.exit()
    else:
        combine_rank(sys.argv[1], sys.argv[2])

