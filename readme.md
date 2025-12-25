--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman 
                                      2025
--------------------------------------------------------------------------------


# PROJECT: Code Rehaul of Too Many Cells Preporcessing Pipline

## 1. OVERVIEW
This program aim is to automate and orgnize the script that was written as part
of the studenst bshara and ___. The pipline produce suitible input for the spectral
clustring algorithm too many cells, which is part of the Immcantation immunulogy 
packages suite.


## 2. PREREQUISITES
Please ensure the following python modules are installed:
- `Pandas`
- `NumPy` 
- `Matplotlib` 
- `SciPy`
- `tqdm`
- `pandarallel`
- `mysql.connector (mysql-connector-python)`


## 3. USAGE GUIDE
1. Congifgure the `congif.json` file (see section 4 - config).
   - Once the custom python modules will be loaded the script will initiate
   the required folders and import the `congif.json` information into `config`
   variable.
2. Run the `main.py` file via cmd

`Note: Run report will be saved in the reports folder`


## 4. CONFIG.JSON CONFIGURATION
The `config.json` is in a json format and it's purpose is to configure the ImmuneDB MySQL connection
and database on which the process will be performed.
before program usage:
- `sql`: Configure the sql connection information.
- `database`: Configure the database name and subject id. 
    - subject_id in the format of "1,2,3,...,n"


## 5. DIRECTORY STRUCTURE
The program uses the following folder structure: 
- `temp_data`: Store the files constructed along the pipeline.
    - `temp_data\{database_subject}`: Dedicated folder for each database and subject.
- `tms_input`: Store the final output, can be used as input for the too many cells algorithm.
    -`tms_input\{database_subject}`: Dedicated folder for each database and subject.
- `reports`: The program will save report for each run to this folder with as `date_database_subject.txt`.


	  
## 6. RESOUCES
- 
