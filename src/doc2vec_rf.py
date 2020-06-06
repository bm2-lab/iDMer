#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: Xiaodi Yang

"""
from __future__ import division
import re, pickle
import os
import numpy as np
import pandas as pd
from os.path import join
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from sklearn.externals.joblib import Parallel
from sklearn.externals.joblib import delayed
import warnings
warnings.filterwarnings('ignore')


def get_documents(seq_list, seq_ids, start, end, k, extract_method):
    """
    :param seq_list:
    :param seq_ids:
    :param start:
    :param end:
    :param k:
    :param extract_method
    Example sequence: MALFFFNNN
    doc2vec parameters: k, extract_method, vector_size, window
    extract_method: [1, 2, 3]
        sequence example: MALFFFNNN
        1) ['MAL', 'ALF', 'LFF', 'FFF', 'FFN', 'FNN', 'NNN']
        2) ['MAL', 'FFF', 'NNN', 'ALF', 'FFN', 'LFF', 'FNN']
        3) ['MAL', 'FFF', 'NNN']
           ['ALF', 'FFN']
           ['LFF', 'FNN']
    vector_size: [32]
    window [3]
    :return:
    """
    from gensim.models.doc2vec import TaggedDocument
    documents = []

    for seq, seq_id in zip(seq_list, seq_ids):
        codes = seq[start: end]
        if extract_method == 1:
            words = [codes[i: i + k] for i in range(len(codes) - (k - 1))]
            documents.append(TaggedDocument(words, tags=[seq_id]))
        elif extract_method == 2:
            words = [codes[j: j + k] for i in range(k) for j in range(i, len(codes) - (k - 1), k)]
            documents.append(TaggedDocument(words, tags=[seq_id]))
        elif extract_method == 3:
            for i in range(k):
                words = [codes[j: j + k] for j in range(i, len(codes) - (k - 1), k)]
                documents.append(TaggedDocument(words, tags=[seq_id + '_%s' % i]))
    return documents


def get_doc2vec_parameters(encode):
    """
    k, extract_method, vector_size, window, epoch
    :param encode:
    :return:
    """
    elements = [int(i) for i in encode.split('-')[2:]]
    try:
        k, extract_method, vector_size, window, epoch = elements
        return k, extract_method, vector_size, window, epoch
    except Exception as e:
        print('Error, the number of hyper-parameter of doc2vec is unequal!!!', elements)
        exit()

def infer_vector(seq_ids, seqs, start, end, encode, model_fname, encoding_fname, complete_infer):
    # get doc2vec parameters
    from gensim.models.doc2vec import Doc2Vec
    k, extract_method, vector_size, window, epoch = get_doc2vec_parameters(encode)
    documents = get_documents(seqs, seq_ids, start, end, k, extract_method)
    model = Doc2Vec.load(model_fname)
    protein_encodings = pickle.load(open(encoding_fname, 'rb'))

    feature_matrix = []
    infernum=0
    unifernum=0
    if complete_infer:
        print('All proteins are inferred!!!')
    for seq_id, document in zip(seq_ids, documents):
        if complete_infer:
            # complete inferring by doc2vec model
            feature_matrix.append(list(model.infer_vector(document[0])))
        else:
            # if seq_id in protein_encodings protein_encodings[seq_id] else inferring by doc2vec model
            if seq_id in protein_encodings:
                feature_matrix.append(protein_encodings[seq_id])
                unifernum+=1
            else:
                print(seq_id, 'is inferred!!!')
                infernum+=1
                feature_matrix.append(list(model.infer_vector(document[0])))
    print('Infered protein number:',infernum)
    print('Uninfered protein number:',unifernum)
    return dict(zip(seq_ids, feature_matrix))


def extract_seq(fasta_fname, min_len):
    seq_ids, seqs = [], []
    for seq_record in SeqIO.parse(fasta_fname, "fasta"):

        bool=re.search(r'[^ACDEFGHIKLMNPQRSTVWY]',(str(seq_record.seq)).upper())
        if bool:continue
        seq=(str(seq_record.seq)).upper()
        if len(seq) > min_len:
            if re.search('\|', seq_record.id):
                seq_ids.append(re.search('\|(\w+)\|', seq_record.id).group(1))
            else:
                seq_ids.append(seq_record.id)
            seqs.append(seq)
        else:
            print(seq_record.id, 'Not encoding')
            continue
    return seq_ids, seqs



def fasta_to_encoding(fasta_fname, encode, start, end, min_len, out_fname, dtype):
    """
    :param fasta_fname:
    :param encode:
    :param start:
    :param end:
    :param out_fname:
    :param dtype:
    :return:
    """
   
    seq_ids, seqs = extract_seq(fasta_fname, min_len)

    # select the encoding method
    feature_matrix, feature_name = encode_select(encode, seqs, start, end,
    seq_ids=seq_ids, out_fname=out_fname)

    if dtype == 'pkl':
        if re.search('\+', encode):
            feature_dict = feature_matrix
        else:
            feature_dict = dict(zip(seq_ids, feature_matrix))
        with open('%s.pkl' % out_fname, 'wb') as f:
            print('fname is ', out_fname)
            pickle.dump(feature_dict, f)
            print(feature_name)
            pickle.dump(feature_name, f)
    elif dtype == 'pd':
        feature_matrix = np.asarray(feature_matrix)
        df = pd.DataFrame(feature_matrix, columns=['f_%s' % (i + 1) for i in range(feature_matrix.shape[1])])
        df.insert(0, 'protein_name', seqs)
        f = open('%s.pd' % out_fname, 'w')
        f.write('# %s\n' % '\t'.join(feature_name))
        df.to_csv(f, sep='\t', index=None)
        f.close()



def get_human_virus_fasta(protein_pair_file,protein_seq_file):
    hvi_candidate=[];pro_seq={}
    r1=open(protein_seq_file)
    for line in r1.readlines():
        line=line.strip()
        pro=line.split()[0]
        seq=line.split()[1]
        pro_seq[pro]=seq
    r1.close()
    
    r2=open(protein_pair_file)
    for line in r2.readlines():
        line=line.strip('\n')
        pro1=line.split('\t')[0]
        pro2=line.split('\t')[1]
        seq1=pro_seq[pro1]
        seq2=pro_seq[pro2]
        if len(seq1)<=30 or len(seq1)>5000:continue
        if len(seq2)<=30 or len(seq2)>5000:continue
        bool1=re.search(r'[^ACDEFGHIKLMNPQRSTVWY]',seq1.upper())
        if bool1:continue
        bool2=re.search(r'[^ACDEFGHIKLMNPQRSTVWY]',seq2.upper())
        if bool2:continue

        hvi_candidate.append('\t'.join([pro1,pro2,seq1,seq2]))

    return(hvi_candidate)

    


def feature_ecode(ecoding,sampleppi_dir,doc2vec_virusshorts,allecoding_proseqfile,ecode,start,end, model_dir, complete_infer,modelname,ecodename):
    import os

    if ecoding=='doc2vec':
    
        seq_ids,seqs=[],[]
        with open(allecoding_proseqfile) as f1:
            for line in f1.readlines():
                line=line.strip('\n')
                pro=line.split('\t')[0]
                seq=line.split('\t')[1]
                seq_ids.append(pro)
                seqs.append(seq)
        f1.close()
        ecoding_dir=sampleppi_dir+'doc2vec/'
        model_fname=model_dir+'human_virus_'+'all'+'-'+ecode+'_'+str(start)+'-'+str(end)+'_'+modelname+'_model.pkl'
        encoding_fname=model_dir+'human_virus_'+'all'+'-'+ecode+'_'+str(start)+'-'+str(end)+'_'+modelname+'.pkl'
        seqids_features=infer_vector(seq_ids, seqs, start, end, ecode, model_fname, encoding_fname, complete_infer)
        for seq_id in seqids_features:
            tmpfeature=[]
            for each in seqids_features[seq_id]:
                tmpfeature.append(str(each))
            seqids_features[seq_id]=tmpfeature


    return(seqids_features)



def get_sample_ecoding(sampleppi_dir,seqids_features,encodinghvpairfile,hvi_candidate):
    import os
    with open(encodinghvpairfile,'w') as w:
        for hvpair in hvi_candidate:
            humanid=hvpair.split('\t')[0]
            virusid=hvpair.split('\t')[1]
            humanseq=hvpair.split('\t')[2]
            virusseq=hvpair.split('\t')[3]
            try:
                w.write(' '.join([humanid,virusid,' '.join(seqids_features[humanid]),' '.join(seqids_features[virusid])])+'\n')
            except:
                print(humanid,virusid,seqids_features)
    
    w.close()

def rf(encodinghvpairfile,rfmodelfile,hvi_candidate_probability,specificity):
    import numpy as np    
    from sklearn.externals import joblib

    data=np.genfromtxt(encodinghvpairfile,dtype=str)
    try:
        ppi_validate=data[:,0:2]
        x=data[:,2:]
        x=x.astype(float)
        num=len(ppi_validate)
    except:
        ppi_validate=data[0:2]
        x=data[2:]
        x=x.astype(float)
        x=x.reshape(1,64)
        num=1


    rfmodel3=joblib.load(rfmodelfile+'3.model')
    prob_predict_y_validate3=rfmodel3.predict_proba(x)
    predictions_validate3=prob_predict_y_validate3[:,1]
        

    with open(hvi_candidate_probability,'w') as w:
        w.write('Pro1ID\tPro2ID\tScore\tInteraction\n')
        for each in range(num):
            if specificity=='0.99' and predictions_validate3[each]>=0.375:
                interaction='yes'
            elif specificity=='0.95' and predictions_validate3[each]>=0.212:
                interaction='yes'
            elif specificity=='0.90' and predictions_validate3[each]>=0.143:
                interaction='yes'
            else:
                interaction='no'
            if num==1:
                w.write(ppi_validate[0]+'\t'+ppi_validate[1]+'\t'+'%.3f'%predictions_validate3+\
                '\t'+interaction+'\n')
            else:
                w.write(ppi_validate[each][0]+'\t'+ppi_validate[each][1]+'\t'+'%.3f'%predictions_validate3[each]+\
                '\t'+interaction+'\n')
       
    w.close()

 

def main():
    import os,sys
    protein_pair_file=sys.argv[1]
    protein_seq_file=sys.argv[2]
    specificity=sys.argv[3]
    Datapath = os.path.dirname(os.path.abspath(__file__))

    sampleppi_dir = sys.argv[4]
    doc2vec_model_dir='{}/doc2vec_model/'.format(Datapath)
    rf_model_dir='{}/rf_model/'.format(Datapath)
    virusshorts=['all']
    allecoding_proseqfile = protein_seq_file
    encodinghvpairfile='{}/Protein_pair_encodings'.format(sampleppi_dir)
    hvi_candidate_probability='{}/PPI_prediction_result.out'.format(sampleppi_dir)
    rfmodelfile=rf_model_dir+'rf_group'
    doc2vec_virusshorts=['all']
    seq_types = ['all-']
    pathogens = ['all']
    seq_ranges = ['0-5000']
    # doc2vec encodes
    doc2vec_encodes = ['doc2vector-%s%s-%s-%s-%s-%s' % (seq_type, k, extract_method, vector_size, window, epoch)
               for seq_type in seq_types for k in [5] for extract_method in [2]
               for vector_size in [32] for window in [3] for epoch in [70]]




    hvi_candidate=get_human_virus_fasta(protein_pair_file,protein_seq_file)


    #doc2vec model
    for ecode in doc2vec_encodes:
        #ecode is 'doc2vector-all-5-2-32-3-70'
        modelname='HVPPI'
        ecodename='HVPPI'

        seqids_features=feature_ecode('doc2vec',sampleppi_dir,['all'],allecoding_proseqfile,ecode,0,5000,doc2vec_model_dir,False,modelname,ecodename)
        get_sample_ecoding(sampleppi_dir,seqids_features,encodinghvpairfile,hvi_candidate)


    rf(encodinghvpairfile,rfmodelfile,hvi_candidate_probability,specificity)


if __name__ == '__main__':
   main()
