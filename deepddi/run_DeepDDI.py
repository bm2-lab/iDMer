import argparse
import os
import shutil
import logging
import time

from deepddi import DeepDDI
from deepddi import preprocessing
from deepddi import result_processing

Datapath = os.path.dirname(os.path.abspath(__file__))
if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-i', '--input_file', required=True, help="Input file")
    parser.add_argument('-p', '--PCA_profile_file', help="PCA profile file")
    
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    
    options = parser.parse_args()
    input_file = options.input_file
    output_dir = options.output_dir
    PCA_profile_file = options.PCA_profile_file

    drug_dir = '{}/data/DrugBank5.0_Approved_drugs/'.format(Datapath)
    pca_model = '{}/data/PCA_tanimoto_model_50.pkl'.format(Datapath)
    trained_weight_model = '{}/data/deepddi_model.h5'.format(Datapath)
    
    multiclass_trained_model = '{}/data/Multiclass_weight.ckpt'.format(Datapath)
    binaryclass_trained_model = '{}/data/Binary_weight.ckpt'.format(Datapath)
    DDI_sentence_information_file = '{}/data/Interaction_information.csv'.format(Datapath)
    
    known_DDI_file = '{}/data/DrugBank_known_ddi.txt'.format(Datapath)
    drug_information_file = '{}/data/Approved_drug_Information.txt'.format(Datapath)
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    ddi_input_file = '%s/tanimoto_PCA50_DDI_Input.csv' % output_dir
    output_file = '%s/DDI_result.txt' % (output_dir)
    ddi_output_file = '%s/Final_DDI_result.txt' % (output_dir)
    annotation_output_file = '%s/Final_annotated_DDI_result.txt' % (output_dir)
    known_drug_similarity_file = '{}/data/drug_similarity.csv'.format(Datapath)
    if PCA_profile_file == None:
        similarity_profile = '%s/similarity_profile.csv' % output_dir
        pca_similarity_profile = '%s/PCA_transformed_similarity_profile.csv' % output_dir
        pca_profile_file = '%s/PCA_transformed_similarity_profile.csv' % output_dir
        print 'calculate structure similarity profile'
        preprocessing.calculate_structure_similarity(drug_dir, input_file, similarity_profile)
        preprocessing.calculate_pca(similarity_profile, pca_similarity_profile, pca_model)
        print 'combine structural similarity profile'
        preprocessing.generate_input_profile(input_file, pca_similarity_profile, ddi_input_file)
    else:
        pca_similarity_profile = PCA_profile_file
        print 'generate input profile'
        preprocessing.generate_input_profile(input_file, pca_similarity_profile, ddi_input_file)
    
    threshold = 0.47
    DeepDDI.predict_DDI(output_dir, output_file, ddi_input_file, trained_weight_model, threshold)    
    result_processing.summarize_prediction_outcome(output_file, ddi_output_file, DDI_sentence_information_file)
    result_processing.annotate_similar_drugs(ddi_output_file, drug_information_file, similarity_profile, known_DDI_file, annotation_output_file, DDI_sentence_information_file, 0.75)
    
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
