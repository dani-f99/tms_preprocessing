# Import required packages
from datetime import datetime
import importlib.metadata
import mysql.connector
import pandas as pd
import unittest
import json
import sys
import os


##############################################################
# Custom function that cheeck if python packages are installed
def check_packages(package_names):
    status = {}

    n_total = len(package_names)
    n_installed = 0
    package_df = pd.DataFrame(index=package_names, columns=["installed","not_installed"], data=0)

    for pkg in package_names:
        try:
            # metadata.version returns the version string if installed
            dist_version = importlib.metadata.version(pkg)
            status[pkg] = f"Installed (v{dist_version})"
            package_df.loc[pkg, "installed"] = 1
            n_installed += 1

        except importlib.metadata.PackageNotFoundError:
            status[pkg] = "Not Installed !!!"
            package_df.loc[pkg, "not_installed"] = 1
    
    print(f"> {n_installed}/{n_total} packages are installed:")

    for i in status:
        print(f"{i} package is {status[i]}")

    return package_df

####################################################################
# DNA codon table - used for the translation from nt to aa sequences
protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
            "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
            "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
            "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
            "TCT" : "s", "CCT" : "P", "ACT" : "T", "GCT" : "A",
            "TCC" : "s", "CCC" : "P", "ACC" : "T", "GCC" : "A",
            "TCA" : "s", "CCA" : "P", "ACA" : "T", "GCA" : "A",
            "TCG" : "s", "CCG" : "P", "ACG" : "T", "GCG" : "A",
            "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
            "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
            "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
            "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
            "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
            "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
            "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
            "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G", 
            "---" : "-"
            }


################################################################################$$##########
# Reading information from json file. Used to extract the parameters from the `config.json`.
def read_json(path:str = "sql_config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open('sql_config.json') as config:
        config_f = json.load(config)

    return config_f


#####################################################
# Creating folder according to the and program scheme
def create_folders(req_folders : list = ["temp_data", "tms_input", "reports"]):
    """
    req_folders : str -> required folders path, if subfolder exsits input '\\' between folders.
    """

    for folder in req_folders:
        if os.path.exists(folder) is False:
            os.mkdir(folder)
            print(f"> folder `{folder}` was created.")

        else:
            print(f"> folder `{folder}` exists, continuing.")


#####################
# sql connector class
class mysql_connector():

    # Getting the connection information from the config
    def __init__(self):
        self.conn_cred = read_json()["sql"]

    # Setting the sql connection
    def setup_conn(self):
        try:
            self.sql_conn = mysql.connector.connect(
                                                    host=self.conn_cred["adress"],
                                                    user=self.conn_cred["username"],
                                                    passwd=self.conn_cred["password"],
                                                    auth_plugin='mysql_native_password',
                                                    )
            print("> Established connection to the MySQL server.")
            
            return self.sql_conn

        except:
            raise Exception("> Failed to establish connection to the MySQL server!")
        
    # Closing the sql connection
    def close_conn(self):
        self.sql_conn.close()
        print("> Connection to the MySQL was closed.")


######################################################
# A small helper to send output to both CMD and a file
class OutputTee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


# ##########################
# Running Pipeline with test
def run_pipeline(test_pipeline,
                 pipeline_name:str = ""
                 ):
    """
    test_pipeline -> the pipeline uninitited unittest pipeline we want to run
    pipeline_name : str -> pipeline name in string format.
    """
    current_time = datetime.now().strftime("%Y-%m-%d-%H-%M")
    db = read_json()["database"]["db_name"]
    reports_path = f"reports\\{db}\\"
    create_folders([reports_path])
    

    # f is your text file -> sys.stdout is the CMD consol
    with open(reports_path+f"{pipeline_name}_[{current_time}]_report_.txt", "w", encoding="utf-8") as f:
        # sys.stdout is the CMD console
        # f is your text file
        dual_stream = OutputTee(sys.stdout, f)
        
        runner = unittest.TextTestRunner(
            stream=dual_stream, 
            verbosity=2, 
            descriptions=True
        )

        # Initialize the runner
        runner = unittest.TextTestRunner(
                stream=dual_stream, 
                verbosity=2, 
                descriptions=True
                )
        
        # unitest  loader object
        loader = unittest.TestLoader()

        # Load tests from the specific class
        suite = loader.loadTestsFromTestCase(test_pipeline) 
        
        # Run with high verbosity for detail
        result = runner.run(suite)
        
        # Custom detailed summary
        print("\n--- PIPELINE EXECUTION SUMMARY ---")
        if result.wasSuccessful():
            print("Final Status: SUCCESS V")
        else:
            print(f"Final Status: FAILED X ({len(result.failures) + len(result.errors)} issues found)")