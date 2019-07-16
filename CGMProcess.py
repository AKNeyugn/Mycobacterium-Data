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
output_first_half = "CGM-strains-A-M.csv"
output_second_half = "CGM-strains-N-Z.csv"

def main():
    start = datetime.datetime.now()
    sys.stdout.write("Start time: " + str(start) + "\n")
    sys.stdout.write("\n")

    input_first = os.path.join(input_folder, first_half)
    input_second = os.path.join(input_folder, second_half)
    output_first = os.path.join(output_folder, output_first_half)
    output_second = os.path.join(output_folder, output_second_half)

    # Create CGM halves
    #process_half(input_first, output_first)
    #process_half(input_second, output_second)

    # Combine CGM halves
    combine_halves(output_first, output_second)

    end = datetime.datetime.now()
    time_taken = end - start
    sys.stdout.write("Time taken: " + str(time_taken.seconds // 60) + " minutes " 
                    + str(time_taken.seconds % 60) + " seconds. \n")
    sys.stdout.write("Script finished! \n")
    return

def process_half(input_name, output_name):
    '''
    Read chemogenomic data file and create CGM output as .csv file

    Args:
        input_name (string): name of chemogenomic data file
        output_name (string): name of CGM .csv file
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
            # Get z_score
            score = row["z_score"]
            cgm_data[strain].append(score)
            sys.stdout.write("Processed %d mutants \r" % (num_cmp))
            sys.stdout.flush()
            num_cmp += 1

    data_df = pd.DataFrame.from_dict(cgm_data)
    indices = dict(zip(list(range(len(cgm_cmps))),cgm_cmps))
    data_df.rename(indices, axis="index", inplace=True)

    data_df.to_csv(output_name)
    sys.stdout.write("Done!\n")
    sys.stdout.write("\n")
    return

def combine_halves(input_one, input_two):
    '''
    Combine the two CGM .csv file halves into one

    Args:
        input_one (string): name of first CGM half
        input_two (string): name of second CGM half
    '''
    sys.stdout.write("Reading CGM halves... \n")
    df_one = pd.read_csv(input_one)
    df_two = pd.read_csv(input_two)
    list_cmps = df_one[df_one.columns[0]]
    df_one.drop(["compound"], axis=1, inplace=True)
    df_two.drop(["compound"], axis=1, inplace=True)

    sys.stdout.write("Combining halves... \n")
    output_name = os.path.join(output_folder, "Full_CGM.csv")
    output_df = pd.concat([df_one, df_two], axis=1, sort=False)
    indices = dict(zip(list(range(len(list_cmps))),list_cmps))
    output_df.rename(indices, axis="index", inplace=True)
    output_df.to_csv(output_name)
    
    sys.stdout.write("Done!\n")
    sys.stdout.write("\n")
    return


if __name__ == "__main__":
    main()