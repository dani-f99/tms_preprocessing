# helper functions import
from .helpers import protein, read_json, create_folders #, mysql_connector

# Python packages import
from sqlalchemy import create_engine, text
from pandarallel import pandarallel
from collections import Counter
from datetime import datetime
from tqdm import tqdm
import pandas as pd
import unittest
import pickle
import math
import csv
import os
import re


class PipelinePreprocessingTest(unittest.TestCase):
    """
    Pipeline of preprocessing ImmuneDB tables for the creation of 
    too many cells input matrix.
    """

    ##################
    # Class initiation 
    @classmethod
    def setUpClass(cls):
        print("--- Initializing Preprocessing Pipeline Environment ---")

        # Importing config information
        cls.config_sql = read_json()["sql"]
        cls.config_db = read_json()["database"]
        cls.db_name = cls.config_db["db_name"]
        cls.time_start = datetime.now()
        cls.subjects = [int(i) for i in cls.config_db["subject_id"].split(",")]
        

        # Cheecking if database assosiated talbes exsists
        create_folders()
        print(cls.db_name, type(cls.db_name))
        create_folders(req_folders=[os.path.join(i, cls.db_name) for i in ["temp_data", "tms_input", "reports"]])

        # Setting paths
        cls.path_temp = os.path.join("temp_data", cls.db_name)
        cls.trimer_dict_path = os.path.join("source", "tables", "trimersDict.csv")


    ################################################################################################
    # Cheecking if the the step already processed the data for this step
    # First step of the preprocessing pipeline, extracting data from the sql server For each subject
    def test_01_translate_convert_extract(self):

        # Itirating over the subjects (list from `json.config`)
        for subject_i in self.subjects:
            toFile = os.path.join(self.path_temp, "1_{}_{}.csv".format(self.db_name, subject_i))
            toFile1 = os.path.join(self.path_temp, "1_AA_{}_{}.csv".format(self.db_name, subject_i))
            toFile2 = os.path.join(self.path_temp, "1_{}_{}_seqK.csv".format(self.db_name, subject_i))

            bool_files = [os.path.exists(i) for i in [toFile, toFile1, toFile2]]

            if sum(bool_files) == 3:
                print("> Step No.1 of the preprocessing already done, continuing to step 2.")
            
            else:
                # Setting-up sql connection credentials
                sql_cred = read_json()["sql"]
                adress, port, username, password = sql_cred.values()
                db_name = read_json()["database"]["db_name"]
                
                # connecting to the MySQL database and defining qry
                engine = create_engine(f"mysql+pymysql://{username}:{password}@{adress}:{port}/{db_name}")
                cmd = f"""SELECT seq.*, coll.* FROM {self.db_name}.sequences 
                                                                           AS seq INNER JOIN {self.db_name}.sequence_collapse 
                                                                           AS coll ON seq.ai=coll.seq_ai WHERE seq.subject_id={int(subject_i)} 
                            AND seq.functional=1 AND coll.instances_in_subject !=0 
                            AND coll.copy_number_in_subject > 1
                            AND seq.deletions is null 
                            AND seq.insertions is null"""

                removed = os.path.join("temp_data", self.db_name, f"{self.db_name}_rem.csv")

                # Getting the data from the sql server
                with engine.connect() as conn:
                    result = conn.execute(text(cmd))
                    seq = [dict(row) for row in result.mappings()]

                #Command for getting the sequences translated:
                print("Translating starts ....")     
                with open(toFile, 'w',newline='') as new_file:
                    csv_writer = csv.writer(new_file)
                    csv_writer.writerow(['seq_id','sequence','TranslatedSeq','TranslatedGermline','ai','subject_id','clone_id','sample_id'])
                    for line in tqdm(seq):
                            dna=""
                            germ=""
                            protein_sequence=""
                            #fix the germline to match the cdr with N's
                            germ=line['germline']
                            cdr3Length=line['cdr3_num_nts']
                            postCDR=line['post_cdr3_length']
                            x=int(cdr3Length)+int(postCDR)
                            replaced=germ[-x:]
                            replaced=replaced.replace('-','N',x) 
                            germ=germ.replace(germ[-x:],replaced)
                            # -------------- To NNNNNNNNNNNNNN #
                            dna=line['sequence']
                            seqID=line['seq_id']
                            ai =line['ai']
                            cloneID=line['clone_id']
                            sampleID=line['sample_id']
                            subjectID=line['subject_id']
                            # Generate protein sequence
                            for i in range(0, len(dna)-(len(dna)%3), 3):
                                if dna[i] == "N" and dna[i+1] == "N" and dna[i+2] == "N":
                                    protein_sequence += "x"
                                elif dna[i] == "N" or dna[i+1] == "N" or dna[i+2] == "N" and dna[i] != "-" and dna[i+1] != "-" and dna[i+2] != "-":
                                    protein_sequence += "x"
                                elif dna[i] == "-" and dna[i+1] == "-" and dna[i+2] == "-":
                                    protein_sequence += protein[dna[i:i+3]]
                                elif dna[i] == "-" or dna[i+1] == "-" or dna[i+2] == "-":
                                    i=i+0;
                                else:
                                    protein_sequence += protein[dna[i:i+3]]
                            germProtein_sequence=""
                            for i in range(0, len(germ)-(len(germ)%3), 3):
                                if germ[i] == "N" and germ[i+1] == "N" and germ[i+2] == "N":
                                    germProtein_sequence += "x"
                                elif germ[i] == "N" or germ[i+1] == "N" or germ[i+2] == "N" and germ[i] != "-" and germ[i+1] != "-" and germ[i+2] != "-":
                                    germProtein_sequence += "x"
                                elif germ[i] == "-" and germ[i+1] == "-" and germ[i+2] == "-":
                                    germProtein_sequence += protein[germ[i:i+3]]
                                elif germ[i] == "-" or germ[i+1] == "-" or germ[i+2] == "-":
                                    i=i+0;
                                else:
                                    germProtein_sequence += protein[germ[i:i+3]]
                            csv_writer.writerow([seqID,dna,protein_sequence,germProtein_sequence,ai,subjectID,cloneID,sampleID])
                    print("Tranlating DONE!")
                
                #Matrix for the mutations
                def build_matrix(rows, cols):
                    matrix = []
                    for r in range(0, rows):
                        matrix.append([0 for c in range(0, cols)])
                    return matrix
                #Mutation function
                def mutatedFunc(seqAA,germAA):
                    global flag
                    flag=0
                    vec=build_matrix(2, len(seqAA))
                    if len(seqAA)!=len(germAA):
                        csv_writer1.writerow([seqAA,germAA])
                        flag=1
                    else:
                        for i in range(0,len(seqAA),1):
                            vec[0][i]=i+1
                        # print(seqAA[i],germAA[i])
                            if seqAA[i]!=germAA[i] and seqAA[i]!= "x" and seqAA[i]!="-" and germAA[i]!= "x" and germAA[i]!="-" and seqAA[i]!="*" and germAA[i]!="*":
                                vec[1][i]=1
                    return vec
                
                print("AA-mutations starts ....")     
                with open(toFile,'r') as csv_file:
                    csv_reader = csv.DictReader(csv_file)
                    with open(toFile1, 'w',newline='') as new_file ,open(removed, 'w',newline='') as nfile:
                        csv_writer = csv.writer(new_file)
                        csv_writer.writerow(['ai','sequence','seq_id','translatedSeq','translatedGerm','vector','subject_id','clone_id','sample_id'])
                        csv_writer1 = csv.writer(nfile)
                        csv_writer1.writerow(['translatedSeq','translatedGerm'])

                        for line in (csv_reader):
                            seq=(line['sequence'])
                            seqID=(line['seq_id'])
                            ai =(line['ai'])
                            seqAA=(line['TranslatedSeq'])
                            germAA=(line['TranslatedGermline'])
                            cloneID=(line['clone_id'])
                            sampleID=(line['sample_id'])
                            subjectID=(line['subject_id'])
                            vec1=mutatedFunc(seqAA, germAA)
                            vector=[]
                            if flag != 1: 
                                for i in range(len(vec1[0])):
                                    if vec1[1][i]==1:
                                        vector.append(i+1)
                                csv_writer.writerow([ai,seq,seqID,seqAA,germAA,vector,subjectID,cloneID,sampleID])
                    print("AA-mutations DONE!")     
                
                def kmersFunc(AA,k):
                    global start
                    start=0
                    p=0
                    kmer=""
                    x=1
                    s=1
                    while(x!=20):
                        if AA[-s] != "-":
                            x+=1
                            s+=1
                        else:
                            s+=1
                
                    for i in range(0,len(AA),1):
                        if AA[i]=="x" or AA[i]=="-":
                            i+=0
                            start+=1
                        else:
                            p1=i
                            for q in range(i,(len(AA)-s)+1,1):
                                for j in range(q,len(AA),1):
                                    if AA[j]=="-" and kmer=="":
                                        j=q+1
                                        p1=j
                                        break
                                    if AA[j]=="-":
                                        j+=0
                                    else:
                                        p+=1
                                        kmer+=AA[j]
                                        if p==k:
                                            p=0
                                            p2=j
                                            pos=(p1+1,p2+1)
                                            #print(AA[p1:p2+1])
                                            #print(kmer)
                                            csv_writer1.writerow([kmer,pos,seqID,ai,subjectID,cloneID,sampleID])
                                            j=q+1
                                            p1=j
                                            kmer=""
                                            break
                            break   
                            
                
                i=0
                print("Kmers extraction starts ....")     
                with open(toFile1,'r') as csv_file:
                    csv_reader = csv.DictReader(csv_file)
                    with open(toFile2, 'w',newline='') as new_file1:
                        csv_writer1 = csv.writer(new_file1)
                        csv_writer1.writerow(['k-mer','position','seq_id','ai','subject_id','clone_id','sample_id'])
                        for line in csv_reader:
                            KmerS=(line['translatedSeq'])
                            seqID=(line['seq_id'])
                            ai =(line['ai'])
                            cloneID=(line['clone_id'])
                            sampleID=(line['sample_id'])
                            subjectID=(line['subject_id'])
                            #Function for the k-mers!
                            kmersFunc(KmerS,20)
                        print("Kmers extraction DONE!")
    
    
    #######################################
    # 2nd step of the preproccsing pipeline
    def test_02_find_unique_by_score(self):
        col_list = ["k-mer","ai","clone_id"]
        DB = self.db_name #database name
        dir_path = self.path_temp #temp file path


        for subject_id in self.subjects:
            output_file = os.path.join(dir_path, "2_{}_{}_byScore.csv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.2 of the preprocessing already done, continuing to step 3.")

            else:
                print("ID = {}".format(subject_id))
                kmers=pd.read_csv(os.path.join(dir_path, "1_{}_{}_seqK.csv".format(DB, subject_id)), usecols=col_list)
                
                kmers = kmers.rename(columns = {'Unique-SeqKmer': 'Kmers'}, inplace = False)
                kmers = kmers.rename(columns = {'k-mer': 'Kmers'}, inplace = False)
                
                
                print(len(kmers))
                kmers['id']=kmers.index
                
                l=kmers.values.tolist()
                
                d = {}
                for i in tqdm(l):
                    d[i[0]] = []
                for j in tqdm(l):
                    if (j[3] not in d[j[0]]):
                        d[j[0]].append(j[3])
                
                
                new_df=pd.DataFrame(kmers['Kmers'])
                l=new_df.values.tolist()
                flat_list = []
                for sublist in l:
                    for item in sublist:
                        flat_list.append(item)
                l = flat_list
                def count_uniqe(lst):
                    new_vals = Counter(l).most_common()
                    new_vals = new_vals[::1] #this sorts the list in scending order
                    return new_vals
                
                new_list=count_uniqe(l)
                
                df = pd.DataFrame(new_list, columns =['kmer', 'score'])

                
                def inx(kmer):
                    return d[kmer][0]
                df['id']=df.kmer.apply(inx)
                
                df.to_csv(output_file, index=False)


    ########################################
    # 3rd step of the preprocessing pipeline
    def test_03_reduce_using_varance_map(self):
        pandarallel.initialize()
        DB = self.db_name #database name
        dir_path = self.path_temp #temp file path
        
        for subject_id in self.subjects:
            input_file = os.path.join(dir_path, "2_{}_{}_byScore.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "3_{}_{}_VarRemain.csv".format(DB,subject_id))

            if os.path.exists(output_file):
                print("> Step No.3 of the preprocessing already done, continuing to step 4.")

            else:
                print("ID = {}".format(subject_id))

                ########################
                col_list = ["kmer","id"]
                col_list1 = ["kmer"]
                kmers0=pd.read_csv(input_file, usecols=col_list)
                kmers1=pd.read_csv(input_file, usecols=col_list1)
                kmers=kmers1.copy()
                kmers1['tmpid']=kmers1.index
                kmers1['var']=0
                
                print(len(kmers0))
                kmers0.head(5)
                
                ########################
                kmers1['id']=kmers0['id']
                
                selected_columns = kmers1[["tmpid","var"]]
                new_df = selected_columns.copy()
                
                all_list = new_df.values.tolist()
                
                All = {}
                st = {}
                for l in tqdm(all_list):
                    All[l[0]] = l[1]
                    st[l[0]] = l[1]
                
                #All
                ########################
                def diff_letters(a,b):
                    cnt=0
                    for i in range(len(a)):
                        if a[i] != b[i] and a[i]!='x' and b[i]!='x':
                            cnt+=1
                            if cnt>1:
                                break
                    return cnt
                
                def wrap_f(x):
                    def fiten(row):
                        row = row[:10]
                        return row
                    return fiten(x['kmer'])
                
                def wrap_l(x):
                    def laten(row):
                        row = row[10:]
                        return row
                    return laten(x['kmer'])
                
                ########################
                kmers['first']=kmers.parallel_apply(wrap_f,axis=1)
                kmers['last']=kmers.parallel_apply(wrap_l,axis=1) 
                kmers['kmer']=kmers.index
                kmers
                
                ########################
                lastlist = kmers.values.tolist()
                
                removed = {}
                for i in tqdm(lastlist):
                    removed[i[0]] = []
                
                last = {}
                for l in tqdm(lastlist):
                    last[l[0]] = l[2]


                first = {}
                for l in tqdm(lastlist):
                    first[l[0]] = l[1]
                

                firstd = {}
                for i in tqdm(lastlist):
                    firstd[i[1]] = []
                
                for j in tqdm(lastlist):
                    firstd[j[1]].append(j[0])


                lastd = {}
                for i in tqdm(lastlist):
                    lastd[i[2]] = []
                
                for j in tqdm(lastlist):
                    lastd[j[2]].append(j[0])

                ########################
                def mapf():
                    r=[]
                    for l in tqdm(firstd.values()):
                        if len(l)>1:
                            for i in range(len(l)):
                                if st[l[i]]!=0:
                                    continue
                                for j in range(i + 1, len(l)):
                                    if st[l[j]]==0:
                                        if diff_letters(last[l[i]],last[l[j]])<=1: 
                                            All[l[i]]+=1
                                            st[l[j]]+=1
                                            removed[l[i]].append(l[j])
                                            r.append(l[j])
                                
                #second comparision
                    for e in tqdm(lastd.values()):
                        if len(e)>1:
                            for x in range(len(e)):
                                if st[e[x]]!=0:
                                    continue
                                for y in range(x + 1, len(e)):
                                    if st[e[y]]==0:
                                        if diff_letters(first[e[x]],first[e[y]])<=1: 
                                            All[e[x]]+=1
                                            st[e[y]]+=1
                                            removed[e[x]].append(e[y])
                                            r.append(e[y])
                    return r
                
                ########################
                trash=mapf()
                ss=list(set(trash))
                print("number of kmers to be deleted: \n{0:d}".format(len(ss)))
                
                ########################
                var=pd.DataFrame.from_dict(All, orient='index',columns=['variance'])
                final = kmers1.drop(columns=['tmpid'])
                final['var'] =var['variance']
                final=final.sort_values(by=['var'],ascending=False)
                final
                
                ########################
                final = final.drop(trash)
                print("the final number of kmers: \n{0:d}".format(len(final)))
                final.head()

                final.to_csv(output_file, index=False)


    #########################################
    # 4th step of the preprocessing  pipeline
    def test_04_siliding_window(self):
        DB = self.db_name #database name
        dir_path = self.path_temp #temp file path
        
        for subject_id in self.subjects:
            input_file = os.path.join(dir_path, "3_{}_{}_VarRemain.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "4_{}_{}_slidingwindow_Var.csv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.4 of the preprocessing already done, continuing to step 5.")

            else:
                df = pd.read_csv(input_file, usecols=["kmer"])
                
                #####
                result=[]
                index_tracker = 0
                threshold = 0
                total_rows = len(df)
                
                #####
                for kmer in df['kmer']:
                    if index_tracker > threshold + 1000000:
                        print("On row {}/{}".format(index_tracker, total_rows))
                        threshold = index_tracker
                    index_tracker += 1
                
                    SlidWindowStr = [kmer[i:i+3] for i in range(len(kmer)-2)]
                    result.append(SlidWindowStr)
                
                df["SlidingWindow"]=result    
                
                #####
                df.to_csv(output_file, index=False)


    ########################################
    # 5th step of the preprocessing pipeline
    def test_05_trimers_filter(self):
        DB = self.db_name #database name
        dir_path = self.path_temp #temp file path
        
        for subject_id in self.subjects:
            input_file = os.path.join(dir_path, "4_{}_{}_slidingwindow_Var.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "5_{}_{}_filtered_trimers_VarRemain.csv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.5 of the preprocessing already done, continuing to step 6.")

            else:
                print("ID = {}".format(subject_id))
    
                # SlidingWindow kmers after Variance PATH
                df = pd.read_csv(input_file)
                
                # Trimers dictionary PATH
                vocab = pd.read_csv(self.trimer_dict_path, index_col=False)
                
                Trimer_dict = pd.Series(vocab.index,index=vocab.trimer).to_dict()
                
                trimers_of_interest = set(vocab['trimer'].tolist())
                
                flag=0
                required_trimers=set()
                not_required_trimers=set()
                for index, row in tqdm(df.iterrows()):
                    kmer_list = row['SlidingWindow']
                    kmer_list = re.sub(r"[^A-Za-z0-9(),]", "", kmer_list)
                    kmer_list = re.sub(r"[^A-Za-z0-9()]", " ", kmer_list)
                    kmer_list = list(kmer_list.split(" "))
                    for kmer in kmer_list:
                        if kmer in trimers_of_interest:
                            required_trimers.add(kmer)
                            flag=1
                        else:
                            not_required_trimers.add(kmer)
                    if flag==0:
                        print(index)
                    flag=0
                
                # required_trimers_sorted
                required_trimers_sorted=[]
                for trimer in vocab['trimer']:
                    if trimer in required_trimers:
                        required_trimers_sorted.append(trimer)
                
                len(required_trimers_sorted)
                
                result=pd.DataFrame(data=required_trimers_sorted, columns=["trimer"])
                
                # Filtered trimers save PATH
                result.to_csv(output_file)


    ########################################
    # 6th step in the preprocessing pipeline
    def test_06_sliding_window_filter(self):
        DB = self.db_name #database name
        dir_path = self.path_temp #temp file path
        
        for subject_id in self.subjects:
            input_df = os.path.join(dir_path, "4_{}_{}_slidingwindow_Var.csv".format(DB, subject_id))
            input_vocab = os.path.join(dir_path, "5_{}_{}_filtered_trimers_VarRemain.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "6_{}_{}_svar_SlidingWindow_filter.csv".format(DB, subject_id))

            if os.path.exists(output_file):
                print("> Step No.6 of the preprocessing already done, continuing to step 7.")

            else:
                # SlidingWindow kmers after variance PATH
                df = pd.read_csv(input_df, usecols=["SlidingWindow"])
                
                vocab = pd.read_csv(input_vocab, index_col=False)
                vocab.drop('Unnamed: 0',axis='columns', inplace=True)
                
                li=set(vocab['trimer'])
                
                reqiured_kmers=set()
                for index, row in tqdm(df.iterrows()):
                    kmer_list = row['SlidingWindow']
                    kmer_list = re.sub(r"[^A-Za-z0-9(),]", "", kmer_list)
                    kmer_list = re.sub(r"[^A-Za-z0-9()]", " ", kmer_list)
                    kmer_list = set(kmer_list.split(" "))
                    for kmer in kmer_list:
                        if kmer in li:
                            reqiured_kmers.add(index)            
                reqiured_kmers

                
                new_list=[]
                for kmer in reqiured_kmers:
                    new_list.append(df.iloc[kmer])
                
                df_new=pd.DataFrame(data=new_list)
                df_new
                
                df_new.loc[df_new.index==349050]
                
                counter=0
                for i in df_new.index:
                    if i != counter:
                        print(counter)
                        counter+=1
                    counter+=1
                
                # Filtered sliding window kmers PATH
                df_new.to_csv(output_file, index=False)


    ###################################################
    # 7th (and last) step of the preprocessing pipeline
    def test_07_finding_trimers_weights(self):
        DB = self.db_name #database name
        dir_path = self.path_temp #temp file path
        
        for subject_id in self.subjects:
            input_df = os.path.join(dir_path, "6_{}_{}_svar_SlidingWindow_filter.csv".format(DB, subject_id))
            output_file = os.path.join(dir_path, "7_{}_{}_VarRemain_trimer_weights.p".format(DB, subject_id))

            # Trimers dictionary PATH
            #Trimers Table    
            TriMers = pd.read_csv(self.trimer_dict_path)

            #Drop Index Column
            TriMers.drop('index',axis='columns', inplace=True)
            # Using DataFrame.insert() to add a column
            TriMers["Places"] = ""
            TriMers

            ### Convert trimer lookup into a dictionary
            trimer_dict = {}
            for index, row in TriMers.iterrows():
                trimer = row['trimer']
                trimer_dict[trimer] = index

            trimer_dict

            if os.path.exists(output_file):
                print("> Step No.7 of the preprocessing already done, run pipeline 2 (matrix creation)")

            else:
                print("ID = {}".format(subject_id))
                
                # Filtered slidingwindow kmers PATH
                #Read CSV File
                df = pd.read_csv(input_df)
                
                ### Find trimer Weights
                #### Make dictionary to save results
                result = {}
                for trimer in TriMers['trimer']:
                    result[trimer] = 0
                
                trimers_of_interest = set(TriMers['trimer'].tolist())
                for index, row in tqdm(df.iterrows()):
                    kmer_list = row['SlidingWindow']
                    kmer_list=eval(kmer_list)
                    #print(type(kmer))
                    for kmer in kmer_list:
                        if kmer in trimers_of_interest:
                            result[kmer]+=1
                
                ## Convert to weight   (IDF Equation)
                total_seqs=len(df) 
                temp={}
                for key,value in result.items():
                    if value!=0:
                        idf=math.log10(total_seqs/value)
                        temp[key]=idf
                result=temp
                
                # Trimers weights Save PATH
                pickle.dump(result, open(output_file, "wb"))