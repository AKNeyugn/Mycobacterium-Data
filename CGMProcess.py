#!/usr/bin/env python

""" Read both chemogenomic data files (A-M and N-Z) and create CGM .csv file:
    If run with *filter_run* argument equal to "true", cell values in CGM will be binary 
    (1 if input values pass z_score AND p_value cutoffs, 0 otherwise); 
    Otherwise, will create 2 CGM files: one where cell values are z_scores and another where cell values are p_values. 

    Author: Roy Nguyen
    Last edited: Aug 2, 2019
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
column_names = ["z_score", "p_value"]

def main():
    start = datetime.datetime.now()
    sys.stdout.write("Start time: " + str(start) + "\n")
    sys.stdout.write("\n")

    input_first = os.path.join(input_folder, first_half)
    input_second = os.path.join(input_folder, second_half)

    filter_run = sys.argv[1]

    if filter_run.upper() == "TRUE":
        z_score_cutoff = sys.argv[2]
        p_value_cutoff = sys.argv[3]

        # Create filtered matrix
        output_name = output_full + "_" + str(z_score_cutoff) + "_" + str(p_value_cutoff) + ".csv"
        # Create DataFrame for CGM halves
        df_first = process_half_filter(input_first, z_score_cutoff, p_value_cutoff)
        df_second = process_half_filter(input_second, z_score_cutoff, p_value_cutoff)
        # Combine CGM halves
        combine_halves(df_first, df_second, output_name)
    else:
        # Create Z-score matrix
        output_name = output_full + "_Z_score.csv"
        # Create DataFrame for CGM halves
        df_first = process_half(input_first, column_names[0])
        df_second = process_half(input_second, column_names[0])
        # Combine CGM halves
        combine_halves(df_first, df_second, output_name)

        # Create Z-score matrix
        output_name = output_full + "_P_value.csv"        
        # Create DataFrame for CGM halves
        df_first = process_half(input_first, column_names[1])
        df_second = process_half(input_second, column_names[1])
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
    Read chemogenomic data file and return DataFrame with desired cell values

    Args:
        input_name (string): name of chemogenomic data file
        cell_value (string): desired column for cell values

    Return:
        (DataFrame): DataFrame with desired cell values
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

def process_half_filter(input_name, z_score_cutoff, p_value_cutoff):
    '''
    Read chemogenomic data file and create DataFrame with binary cell values;
    1 if input values pass z_score AND p_value cutoffs, 0 otherwise 

    Args:
        input_name (string): name of chemogenomic data file
        z_score_cutoff (float): z_score cutoff
        p_value_cutoff (float): p_value cutoff
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
            z_score = row[column_names[0]]
            p_value = row[column_names[1]]
            if float(z_score) <= float(z_score_cutoff) and float(p_value) <= float(p_value_cutoff):
                score = "1"
            else:
                score = "0"
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