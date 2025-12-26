# Imports
from source.helpers import run_pipeline, check_packages
from source.pipeline1_preprocessing import PipelinePreprocessingTest
from source.pipeline2_matrixmaker import PipelineMatrixMakerTest
import pandas as pd

# Cheeking that the required python packages are installed
installed_packages = check_packages(["pandas", "numpy", "scipy", "tqdm", "pandarallel", "mysql-connector-python"])
missing = installed_packages.loc[installed_packages.not_installed == 1, "not_installed"].index.to_list()

# If some missing packages exists -> raise error
if len(missing) > 0:
    raise Exception(f"\nplease install required packages: {missing}")

# Execute the pipeline
if __name__ == "__main__":
    pipeline1_result = run_pipeline(PipelinePreprocessingTest, pipeline_name="preprocessing_pipeline1")
    pipeline2_results = run_pipeline(PipelineMatrixMakerTest, pipeline_name="matrixmaker_pipeline2")