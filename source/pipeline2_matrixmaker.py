# Custom functions import
from .helpers import read_json, mysql_qry

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
        cls.subjects = [int(i) for i in cls.config_db["subject_id"].split(",")]
        

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
            output_file = os.path.join(dir_path, "{}_{}_matrix.mtx".format(DB, subject_id))

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
                shape_rows = len(Kmers)
                shape_columns = len(Trimers)
                print("The shape of the dense matrix is: {} rows x {} columns".format(shape_rows, shape_columns))
                
                row = np.array(Kmers)
                col = np.array(list(Trimers))
                data = np.array(Weights)
                
                sparse_matrix=csr_matrix((data, (col,row)), shape=(shape_columns,shape_rows))
                
                # Sparse matrix save PATH  
                scipy.io.mmwrite(output_file, sparse_matrix)
    
    
    #############################################################################
    # 2nd step of the tsm input preparation - barcodes creation (k-mers, columns)
    def test_02_barcodes_creator(self):
        DB = self.db_name #database name
        temp_path = self.path_temp # first preprocessing pipeline path
        dir_path = self.path_final #final file path
        # brcodes (kmers) file: 6_{}_{}_filt_slidingwindow_Var.csv
        
        for subject_id in self.subjects:
            input_temp_kmer_path = os.path.join(temp_path, "6_{}_{}_filt_slidingwindow_Var.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "{}_{}_barcodes.tsv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.2 of the matrix creation pipeline already done, continuing to step 3.")
            
            else:
                temp_kmer_df = pd.read_csv(input_temp_kmer_path, index_col=1)
                temp_kmer_df_kmers = temp_kmer_df.index
                
                with open(output_file, 'wt') as out_file:
                    tsv_writer=csv.writer(out_file,delimiter='\t')
                    for i in range(len(temp_kmer_df_kmers)):
                        tsv_writer.writerow([temp_kmer_df_kmers[i]])

    
    #############################################################################
    # 3rd step of the tsm input preparation - genes list creation (trimers, rows)
    # genes (trimers) file: "5_{}_{}_filtered_trimers_VarRemain.csv"
    def test_03_feature_creator(self):
        DB = self.db_name #database name
        temp_path = self.path_temp # first preprocessing pipeline path
        dir_path = self.path_final #final file path
        
        for subject_id in self.subjects:
            input_temp_trimer_path = os.path.join(temp_path, "5_{}_{}_filtered_trimers_VarRemain.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "{}_{}_genes.tsv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.3 of the matrix creation pipeline already done, process is done.")
            
            else:
                temp_trimer_df = pd.read_csv(input_temp_trimer_path, index_col=0)
                temp_trimer_list = temp_trimer_df.trimer.values
                
                with open(output_file, 'wt') as out_file:
                    tsv_writer=csv.writer(out_file,delimiter='\t')
                    
                    for i in range(len(temp_trimer_list)):
                        tsv_writer.writerow([temp_trimer_list[i]])

    #########################################################
    # 4th step of the tsm input preparation - labels creation
    # brcodes (kmers) file: 6_{}_{}_filt_slidingwindow_Var.csv
    def test_04_labels_creator(self):
        DB = self.db_name #database name
        temp_path = self.path_temp # first preprocessing pipeline path
        dir_path = self.path_final #final file path

        for subject_id in self.subjects:
            output_file = os.path.join(dir_path, "{}_{}_labels.csv").format(DB, subject_id) # path for the labels.csv

            if os.path.exists(output_file):
                print("> Step No.3 of the matrix creation pipeline already done, process is done.")
            
            else:
                # Importing the metadata table from the MySQL serever
                # Setting up label and mysql query
                metadata_label = self.config_db["metadata_label"]
                metadata_qry = "SELECT * FROM covid_vaccine_new.sample_metadata"
                metadata_df = mysql_qry(qry=metadata_qry)
                
                
                # Verifing that there is `metadata_label` input in the `sql_congif.json` file
                if metadata_label == "":
                    print("> No metadata label selected, aborting creation of `labels.csv`")
                
                else:
                    print(f"metadata label `{metadata_label}` selected, creating `labels.csv`.")

                    # Dictionary that link sample_id to metadata label:
                    filt_df = metadata_df[metadata_df["key"] == metadata_label]
                    labels_dict = {int(i):j for i,j in zip(filt_df.sample_id.values, filt_df.value.values)}

                    # Loading the required datasets in order backtrack to the first metadata contining files
                    """
                    Steps of the backtracking:
                    1. getting the `id` labels from the k-mers used to create the barcodes (f03).
                    2. merging the filterd k-mers (f03) with the k-mer table (f01) by 'id'.
                    """
                
                    # Importing the preprocsssed tables
                    input_f03_kmers = os.path.join(temp_path, "3_{}_{}_VarRemain.csv".format(DB, subject_id)) # kmers files with id column
                    input_f06_filt_kmers = os.path.join(temp_path, "6_{}_{}_filt_slidingwindow_Var.csv".format(DB, subject_id)) # filtred k-mer file (used for barcodes creation)
                    input_01_seqk = os.path.join(temp_path, "1_{}_{}_seqK.csv".format(DB, subject_id)) # initial k-mer file (with metadata)
                    input_barcodes = os.path.join(dir_path, "covid_vaccine_new_7_barcodes.tsv".format(DB, subject_id)) # path of the barcodes.tsv

                   
                    # Joining the tables to get the sample id
                    f01_df = pd.read_csv(input_01_seqk).reset_index(names="id")
                    f03_df = pd.read_csv(input_f03_kmers)
                    f06_df = pd.read_csv(input_f06_filt_kmers, index_col=0)
                    barcodes_df = pd.read_csv(input_barcodes, sep="\t", index_col=None, header=None)
                    barcodes_df.columns = ["kmer"]

                    # Getting the filtred kmers f03
                    f03_df = f03_df[f03_df.index.isin(f06_df.kmer.values)]

                    # Using the labels_dict to map the labels (via sample id)
                    merged_df = pd.merge(left=f01_df, right=f03_df, how="right", on="id")[["kmer","id","sample_id"]]
                    merged_df["label"] = merged_df.sample_id.map(labels_dict)

                    # Adding labels to the k-kmers
                    labels_df = pd.merge(left=barcodes_df, right=merged_df, on="kmer", how="left")[["kmer", "label"]]
                    labels_df.columns = [["item","label"]]
                    labels_df.to_csv(output_file, index=0)

                    print(f"{output_file} labels file created.")
