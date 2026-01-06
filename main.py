# Imports
from source.helpers import run_pipeline, check_packages, create_folders, read_json
from source.pipeline1_preprocessing import PipelinePreprocessingTest
from source.pipeline2_matrixmaker import PipelineMatrixMakerTest
import os


# Cheeking that the required python packages are installed
installed_packages = check_packages(["pandas", "numpy", "scipy", "tqdm", "pandarallel", "mysql-connector-python"])
missing = installed_packages.loc[installed_packages.not_installed == 1, "not_installed"].index.to_list()


# If some missing packages exists -> raise error
if len(missing) > 0:
    raise Exception(f"\nplease install required packages: {missing}")
    

# Create required folders -> the input is the built in argument input, it's written here for clearity
# getting database and subject information from config file `sql_config.json`
db_info = read_json()["database"]
database_name, subjects = db_info["db_name"], db_info["subject_id"].split(",")


# Creating folders
main_folders = ["temp_data", "tms_input", "reports"] # main folders
sub_folders = [os.path.join("tms_input", f"{database_name}-subject{i}") for i in subjects] # sub
create_folders(main_folders + sub_folders) 


# Execute the pipeline
if __name__ == "__main__":
    pipeline1_result = run_pipeline(PipelinePreprocessingTest, pipeline_name="preprocessing_pipeline1")
    pipeline2_results = run_pipeline(PipelineMatrixMakerTest, pipeline_name="matrixmaker_pipeline2")