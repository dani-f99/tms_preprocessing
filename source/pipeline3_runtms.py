from .helpers import read_json, create_folders
import unittest
import os

# Pipeline for too many cells run
class RunTmsTest(unittest.TestCase):
    """
    Pipeline that will run the too many cells spectral clustring analysis
    on the pipeline 2 output (pipeline2_matrixmaker).
    """

    ##################
    # Class initiation 
    @classmethod
    def setUpClass(cls):
        print("--- Initializing Too Many Cells Pipeline Environment ---")
        
        # Importing database and subjects information -> for folder creation
        config = read_json()
        cls.db_name = config["database"]["db_name"]
        cls.db_subjects = config["database"]["subject_id"].split(",")

    # Creating required folders
    def test_01_folders_creation(self):        
        path_mainf = os.path.join("tms_output", self.db_name)
        req_folders = ["tms_output", path_mainf] + [os.path.join(path_mainf,i) for i in self.db_subjects]
        create_folders(req_folders)











































def run_tms_pipeline():
    # 1. Load configuration from sql_config.json
    config_path = Path("sql_config.json")
    if not config_path.exists():
        print(f"Error: {config_path} not found.")
        return

    with open(config_path, "r") as f:
        config = json.load(f)

    db_name = config["database"]["db_name"]
    # Convert subject_id string "1,2,3" into a list ["1", "2", "3"]
    subject_ids = [s.strip() for s in config["database"]["subject_id"].split(",")]

    # 2. Define main paths
    main_dir = Path.cwd()
    input_root = main_dir / "tms_input"
    output_root = main_dir / "tms_output"

    # Create tms_output folder if it doesn't exist
    output_root.mkdir(parents=True, exist_ok=True)

    # 3. Iterate over subjects
    for sub_id in subject_ids:
        folder_name = f"{db_name}-{sub_id}"
        
        # Define specific subject paths
        sub_input_dir = input_root / folder_name
        sub_output_dir = output_root / folder_name
        
        # Create the subject-specific output directory
        sub_output_dir.mkdir(parents=True, exist_ok=True)
        
        labels_file = sub_input_dir / "labels.csv"
        matrix_path = sub_input_dir / "input" # Assuming 'input' matrix is inside the subject folder

        if not labels_file.exists():
            print(f"Warning: Labels file not found for {folder_name}. Skipping...")
            continue

        print(f"Processing: {folder_name}...")

        # 4. Construct and Run the TMS command
        # The command is executed with sub_output_dir as the working directory
        tms_cmd = [
            "too-many-cells", "make-tree",
            "--matrix-path", str(matrix_path),
            "--labels-file", str(labels_file),
            "--draw-collection", "PieRing",
            "--output", "out"
        ]

        try:
            subprocess.run(
                tms_cmd, 
                cwd=str(sub_output_dir), # Run command inside the specific output folder
                check=True
            )
            print(f"Successfully processed {folder_name}")
            
        except subprocess.CalledProcessError as e:
            print(f"Error running TMS for {folder_name}: {e}")

if __name__ == "__main__":
    run_tms_pipeline()