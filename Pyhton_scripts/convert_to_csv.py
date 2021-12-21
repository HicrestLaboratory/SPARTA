# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:16:49 2021

@author: Paolo
"""


import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse

def convert_to_CSV(filename, new_file_name):
    with open(filename, 'r') as infile:
        with open(new_file_name, 'w') as outfile:
            fieldsline = infile.readline();
            
            firstrow = ",".join(fieldsline.split(" "))
            
            outfile.writelines(firstrow);
        
            stop = False;
            count = 0;
            while ( not stop ):
                    line = infile.readline();
                    if not line:
                        stop = True;
                        break
                    newline = ",".join(line.split(" "))
                    outfile.writelines(newline);                    
                    line = infile.readline();
                    count+= 1;
            print("read ", count, " lines");
                    
                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes csv from experiments")
    
    parser.add_argument("--input-file", default="/",
        help="file that contains the already run experiments")
    parser.add_argument("--output-csv", default = "",
        help="file of the new csv to output")
    

    args = parser.parse_args()

    input_file = args.input_file;
    output_csv = args.output_csv;
    if output_csv == "":
        output_csv = input_file.split(".")[0] + ".csv"
                    
convert_to_CSV(input_file, output_csv);
results_df = pd.read_csv(output_csv); #check that results can be read

