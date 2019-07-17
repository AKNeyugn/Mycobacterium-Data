# Mycobacterium-Data

**Paper:** Johnson et al. 2019. "Large-Scale Chemical–Genetics Yields New *M. tuberculosis* Inhibitor Classes." *Nature* 571 (July): 72-78.
https://doi.org/10.1038/s41586-019-1315-z

# CompoundStructure.py

For each compound in given library, convert SMILES string into 3D structure and output .pdb file.

**Module requirements:**

- OpenBabel (http://openbabel.org/wiki/Main_Page)
- pandas (pip install pandas)

**How to run:** *python CompoundStructure.py library_SMILES_file single_file*
- library_SMILES_file: path from current working directory to SMILES file of library to process (eg. Data-Files\compounds.csv)
- single_file: "true" if want output to be single file containing all compounds info for each library

**Outputs:**
- *Compound-3D-Structure* folder containing seperate .pdb files of each compounds in given library
- if running with single_file = "true", *Compound-3D-Structure* folder will also contain .pdb file containing all compounds info for each library

**Note:**
Expected compounds with failed 3D structure conversion due to missing SMILES:
- U14272896
- U20739209
- U21835532
- U32963902
- U84341231
- nan

The files for these compounds are manually moved to seperate *nan-SMILES* folder

**Note:**
.pdb file containing all compounds info is not uploaded to repository due to its large size

# CGMProcess.py

Read both chemogenomic data files (A-M and N-Z) and create 2 seperate chemogenomic matrix (CGM) outputs as .csv files, then combine the halves into one full CGM .csv file.

**Module requirements:**

- pandas (pip install pandas)

**How to run:** *python CGMProcess.py*

**Outputs:**
- *CGM* folder containing 3 .csv files: 2 incomplete CGM files (A-M and N-Z) and 1 full CGM file that combines the previous 2

**Note:**
CGM output is a compound vs gene matrix, where x-axis = gene and y-axis = compound

**Note:**

*A-M (input):*

- Total # rows: 4636478
- Unique compound rows with non n/a compound_stem: 47272
- Unique strain rows: 98

*N-Z (input):*

- Total # rows: 2696727
- Unique compound rows with non n/a compound_stem: 47272
- Unique strain rows: 57

*Full (output):*

- Total # rows (compounds): 47272
- Total # columns (strains): 155
- Unique compound rows with non n/a compound_stem: 47272

**Note:**
Input files, as well as output files, are not uploaded to repository due to their large size
