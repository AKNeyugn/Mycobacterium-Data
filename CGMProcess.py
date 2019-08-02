#!/usr/bin/env python

""" Read both chemogenomic data files (A-M and N-Z) 
    and create 2 seperate chemogenomic matrix (CGM) outputs as .csv files,
    then combine the halves into one full CGM .csv file

    Author: Roy Nguyen
    Last edited: July 16, 2019
"""


import sys
import os
import datetime
import time
import json
import pandas as pd

input_folder = "Data-Files"
output_folder = "CGM"
first_half = "2019-05-15_split-lfc-p-values-strains-A-M.csv"
second_half = "2019-05-15_split-lfc-p-values-strains-N-Z.csv"
output_full = "Full_CGM"

def main():
    start = datetime.datetime.now()
    sys.stdout.write("Start time: " + str(start) + "\n")
    sys.stdout.write("\n")

    input_first = os.path.join(input_folder, first_half)
    input_second = os.path.join(input_folder, second_half)

    filter_run = sys.argv[1]

    if filter_run.upper() == "FALSE":
        # Create Z-score matrix
        output_name = output_full + "_Z_score.csv"
        # Create DataFrame for CGM halves
        df_first = process_half(input_first, "z_score")
        df_second = process_half(input_second, "z_score")
        # Combine CGM halves
        combine_halves(df_first, df_second, output_name)

        # Create Z-score matrix
        output_name = output_full + "_P_value.csv"        
        # Create DataFrame for CGM halves
        df_first = process_half(input_first, "p_value")
        df_second = process_half(input_second, "p_value")
        # Combine CGM halves
        combine_halves(df_first, df_second, output_name)

    end = datetime.datetime.now()
    time_taken = end - start
    sys.stdout.write("Time taken: " + str(time_taken.seconds // 60) + " minutes " 
                    + str(time_taken.seconds % 60) + " seconds. \n")
    sys.stdout.write("Script finished! \n")
    return

def process_half(input_name, cell_value):
    '''
    Read chemogenomic data file and create CGM output as .csv file

    Args:
        input_name (string): name of chemogenomic data file
        cell_value (string): desired column for cell values
    '''
    sys.stdout.write("Processing " + input_name + "...\n")
    df = pd.read_csv(input_name)
    cgm_data = {}
    cgm_cmps = []

    for index, row in df.iterrows():
        strain = row["strain"]
        if strain not in cgm_data.keys():
            # break to shorten output data during testing
            #if len(cgm_data.keys()) == 2:
                #break
            cgm_data[strain] = []
            num_cmp = 1
            sys.stdout.write("\n")
            sys.stdout.write("Processing " + str(strain) + "...\n")

        cmp_stem = row["compound_stem"]
        if str(cmp_stem) != "nan":
            # Get compound name
            compound = row["compound"]
            conct = row ["concentration"]
            cmp_name = compound + "-" + str(conct)
            if cmp_name not in cgm_cmps:
                cgm_cmps.append(cmp_name)
            # Get cell values
            score = row[cell_value]
            cgm_data[strain].append(score)
            sys.stdout.write("Processed %d mutants \r" % (num_cmp))
            sys.stdout.flush()
            num_cmp += 1

    data_df = pd.DataFrame.from_dict(cgm_data)
    indices = dict(zip(list(range(len(cgm_cmps))),cgm_cmps))
    data_df.rename(indices, axis="index", inplace=True)

    sys.stdout.write("Done!\n")
    sys.stdout.write("\n")
    return data_df

def combine_halves(df_one, df_two, output):
    '''
    Combine the two CGM .csv file halves into one

    Args:
        df_one (string): DataFrame of first CGM half
        df_two (string): DataFrame of second CGM half
        output (string): name of combined CGM
    '''
    list_cmps = df_one[df_one.columns[0]]
    df_one.drop(["compound"], axis=1, inplace=True)
    df_two.drop(["compound"], axis=1, inplace=True)

    sys.stdout.write("Combining halves... \n")
    output_name = os.path.join(output_folder, output)
    output_df = pd.concat([df_one, df_two], axis=1, sort=False)
    indices = dict(zip(list(range(len(list_cmps))),list_cmps))
    output_df.rename(indices, axis="index", inplace=True)
    output_df.to_csv(output_name)
    
    sys.stdout.write("Done!\n")
    sys.stdout.write("\n")
    return


if __name__ == "__main__":
    main()