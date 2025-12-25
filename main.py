# Imports
from source.helpers import run_pipeline
from source.pipeline1_preprocessing import PipelinePreprocessingTest
from source.pipeline2_matrixmaker import PipelineMatrixMakerTest

# Execute the pipeline
if __name__ == "__main__":
    pipeline1_result = run_pipeline(PipelinePreprocessingTest, pipeline_name="preprocessing_pipeline1")
    pipeline2_results = run_pipeline(PipelineMatrixMakerTest, pipeline_name="matrixmaker_pipeline2")