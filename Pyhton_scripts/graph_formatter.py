# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:59:00 2020

@author: Paolo
"""

import os

os.chdir("../data/")
documents = os.listdir()

selected_idx = 0;
while True:
    print("Select a file:")
    for i, doc in enumerate(documents):
        print(i, doc)
    
    selected_idx = input("Select a file to fix")
    if selected_idx == "x":
        exit()
    if not selected_idx.isnumeric():
        print("wrong selection")
        continue;
    selected_idx = int(selected_idx);
    selected_doc = documents[selected_idx]
    
    with open(selected_doc, 'r') as f:
        line = "#"
        print(line[0])
        while not (line[0].isnumeric()):
            line = f.readline();
        i = 0
        while((line[i]).isnumeric() or line[i] == '.'):
            i += 1;
        separator = line[i];
        
    print("separator: ", repr(separator))
    new_separator = input("Insert new separator:")
    if (new_separator == "\n"):
        print("exiting")
        exit()
    
    with open(selected_doc, 'r') as f:
        with open("fixed_" + selected_doc, 'w') as write_f:
            for line in f:
                line = f.readline();
                if (not (line[0]).isnumeric()):
                    continue
                splitline = line.split(separator);
                newline = splitline[0] + new_separator + splitline[1] + "\n";
                write_f.write(newline);
            
    print("new file created")