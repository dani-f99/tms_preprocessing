import json


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
def read_json(path:str = "config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open('config.json') as config:
        config_f = json.load(config)

    return config_f


#####################################################
# Creating folder according to the and program scheme
def create_folders():
    req_folders = ["temp_data", "tms_input", "reports"]

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