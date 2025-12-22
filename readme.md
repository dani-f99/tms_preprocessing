--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman
                                      2025
--------------------------------------------------------------------------------


# PROJECT: Amino Acid Motif Logo Maker

## 1. OVERVIEW
This program aim is to visualize the amino acid usage bais across 
defiened BCR heavy chain variable section. If metadata exists the program
can output multiple sub-plots defined by the metadata in order to comapre 
motifs across different conditions.


## 2. PREREQUISITES
Please ensure the following python modules are installed:
- `Pandas`
- `NumPy` 
- `Matplotlib` 
- `logomaker`


## 3. USAGE GUIDE
1. Congifgure the `congif.json` file (see section 4 - config).
   - Once the custom python modules will be loaded the script will initiate
   the required folders and import the `congif.json` information into `config`
   variable.
1. Place the sequences dataset in the `input` folder (example file is provided). 
2. Follow the steps illustrated in the `motif_logo.ipynb` notebook(example output 
is provided). 
3. Access the output figures via the `output` folder.


## 4. CONFIG.JSON CONFIGURATION
The `congif.json` is in a json format and serve multiple purposes, configure 
before program usage:
- `input_folder`: raw sequences inputs folder (defualt `input`).
- `output_folder`:raw sequences output folder (defualt `output`).
- `group_by`: column on which the data will be grouped by.
- `split_subjects`: Split the sequences per subject.


## 5. DIRECTORY STRUCTURE
The program uses the following folder structure: 
- `input/`  : Input folder (for the sequences data). 
- `output/` : Result - main output folder. 
  - `output/motif_figure` : output folder of motif figures. 
  - `output/motif_data` : output folder of motif processed data. 

		   
## 6. RESOUCES
- Logomaker documentation: https://logomaker.readthedocs.io/en/latest/
