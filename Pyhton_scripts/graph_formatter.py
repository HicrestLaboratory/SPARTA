# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:59:00 2020

@author: Paolo
"""

import os

os.chdir("../data/")
documents = os.listdir()

print("Select a file:")
for i, doc in enumerate(documents):
    print(i, doc)

selected_idx = int(input())
selected_doc = documents[selected_idx]

with open(selected_doc, 'r') as f:
    with open("fixed_" + selected_doc, 'w') as write_f:
        for line in f:
            line = f.readline();
            splitline = line.split();
            newline = splitline[0] + "\t" + splitline[1] + "\n";
            write_f.write(newline);
        