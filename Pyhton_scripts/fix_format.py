# -*- coding: utf-8 -*-
"""
Created on Wed May  5 11:58:38 2021

@author: Paolo
"""

datasetName = "C:/Users/Paolo/OneDrive/Documenti/GitHub/SPARTA/results/reordering_vs_saad"

with open(datasetName + ".txt", 'r') as inf:
    with open(datasetName + "_fixed.txt", 'w') as outf:
        while (True):
            line = inf.readline();
            if not line:
                break;
            if (line[0] != "r"):
                outf.writelines(line);