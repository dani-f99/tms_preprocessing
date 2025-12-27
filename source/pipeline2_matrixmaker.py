# Custom functions import
from .helpers import read_json

# Import of python packages
from scipy.sparse import csr_matrix
from datetime import datetime
from tqdm import tqdm
import pandas as pd
import numpy as np
import unittest
import scipy
import csv
import re
import os


# Matrix creation class
class PipelineMatrixMakerTest(unittest.TestCase):
    """
    Pipeline of matrix creation, this is the second pipeline of the preprocessing
    required for the tms input file.
    """

    ##################
    # Class initiation 
    @classmethod
    def setUpClass(cls):
        print("--- Initializing MatrixMaker Pipeline Environment ---")

        # Importing config information
        cls.config_db = read_json()["database"]
        cls.db_name = cls.config_db["db_name"]
        cls.time_start = datetime.now()
        cls.subjects = list(cls.config_db["subject_id"])
        

        # Setting paths
        
        cls.path_temp = os.path.join("temp_data", cls.db_name)
        cls.path_final = os.path.join("tms_input", cls.db_name)


    ################################
    # 1st step in the tms input data
    def test_01_matrix_builder(self):
        DB = self.db_name #database name
        temp_path = self.path_temp
        dir_path = self.path_final #final file path
        
        for subject_id in self.subjects:
            input_sliding_window = os.path.join(temp_path, "6_{}_{}_svar_SlidingWindow_filter.csv".format(DB, subject_id))
            input_trimers_weights = os.path.join(temp_path, "7_{}_{}_VarRemain_trimer_weights.p".format(DB, subject_id))
            input_vocab = os.path.join(temp_path, "5_{}_{}_filtered_trimers_VarRemain.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "\\{}_{}_matrix".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.1 of the matrix creation pipeline already done, continuing to step 2.")

            else:
                print("id = {}".format(subject_id))
    
                # covidVacc_id Data SlidingWindow
                df = pd.read_csv(input_sliding_window, usecols=["SlidingWindow"])
                
                # Trimers Table   (IDF)
                trimer_weights = pd.read_pickle(input_trimers_weights)
                
                # Vocab (filtered trimers)
                vocab = pd.read_csv(input_vocab, index_col=False)
                vocab.drop('Unnamed: 0',axis='columns', inplace=True)
                #_VarRemain.csv
                
                # Trimers Dict
                Trimer_dict=pd.Series(vocab.index,index=vocab.trimer).to_dict()
                
                # Lists
                trimers_of_interest = set(vocab['trimer'].tolist())
                Trimers=[]
                Kmers=[]
                Weights=[]
                trimers_counter=0
                kmers_counter=0
                length_counter=0
                
                for index, row in tqdm(df.iterrows()):
                    kmer_list = row['SlidingWindow']
                    kmer_list = re.sub(r"[^A-Za-z0-9(),]", "", kmer_list)
                    kmer_list = re.sub(r"[^A-Za-z0-9()]", " ", kmer_list)
                    kmer_list = set(kmer_list.split(" "))
                    for kmer in kmer_list:
                        if kmer in trimers_of_interest:
                            if Trimer_dict[kmer] not in Trimers:
                                trimers_counter+=1
                            Trimers.append(Trimer_dict[kmer])
                            Kmers.append(kmers_counter)
                            Weights.append(trimer_weights[kmer])
                            length_counter+=1
                            
                    kmers_counter+=1
                shape_rows = max(Kmers)+1
                shape_columns = max(Trimers)+1
                print("The shape of the dense matrix is: {} rows x {} columns".format(shape_rows, shape_columns))
                
                row = np.array(Kmers)
                col = np.array(list(Trimers))
                data = np.array(Weights)
                
                sparse_matrix=csr_matrix((data, (col,row)), shape=(shape_columns,shape_rows))
                
                # Sparse matrix save PATH  
                scipy.io.mmwrite(output_file, sparse_matrix)
    
    
    #######################################
    # 2nd step of the tsm input preparation
    def test_02_barcodes_creator(self):
        DB = self.db_name #database name
        temp_path = self.path_temp
        dir_path = self.path_final #final file path
        
        for subject_id in self.subjects:
            input_temp_kmer_path = os.path.join(temp_path, "3_{}_{}_VarRemain.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "\\{}_{}_barcodes.tsv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.2 of the matrix creation pipeline already done, continuing to step 3.")
            
            else:
                temp_kmer_path = input_temp_kmer_path
                temp_kmer = pd.read_csv(temp_kmer_path).shape[0] + 1
                
                with open(output_file, 'wt') as out_file:
                    tsv_writer=csv.writer(out_file,delimiter='\t')
                    for i in range(1,temp_kmer):
                        tsv_writer.writerow([i])

    
    #######################################
    # 3rd step of the tsm input preparation
    def test_03_feature_creator(self):
        DB = self.db_name #database name
        temp_path = self.path_temp
        dir_path = self.path_final #final file path
        
        for subject_id in self.subjects:
            input_temp_trimer_path = os.path.join(temp_path, "5_{}_{}_filtered_trimers_VarRemain.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "\\{}_{}_genes.tsv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.3 of the matrix creation pipeline already done, process is done.")
            
            else:
                temp_trimer = pd.read_csv(input_temp_trimer_path).shape[0]+1
                
                with open(output_file, 'wt') as out_file:
                    tsv_writer=csv.writer(out_file,delimiter='\t')
                    
                    for i in range(1,temp_trimer):
                        tsv_writer.writerow([i])