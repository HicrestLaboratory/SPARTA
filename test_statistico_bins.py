# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 17:04:21 2020

@author: Paolo
"""
import random as rn
import matplotlib.pyplot as plt
import os

saving_path = "images/reorder_experiment";
try:
    os.mkdir(saving_path)
except:
    print("overwriting images");

sparsities = [0.001,0.01,0.03, 0.05,];

plt.figure()
for s in sparsities:
    
    l = 1000
    a = [0 for i in range(l)];
    
    nz = int(l*s);
    
    for i in range(nz):
        a[i] = 1;
        
    rn.shuffle(a);
    
    x =  [1,5,10,20,30,50,64,80,120, 256];
    y = [];
    y_theory = []
    for block_size in x:
        block_n = int(l/block_size);
    
        print("size", block_size, "n", block_n);
        
        
        full_count = 0;
        for bi in range(block_n):
            for i in range(block_size):
                if a[bi*block_size + i] != 0:
                    full_count += 1;
                    break;
        y.append(full_count/block_n);
        
        yt = ((block_n - 1)/block_n)**(l*s)
        yt = 1 - yt;
        yt = s/yt
        y_theory.append(yt);
        
    
    plt.ylim([0,0.15])
    plt.plot(x,y_theory, label = s)     
    plt.ylabel("density in blocks");
    plt.xlabel("size of blocks");
    
    plt.legend(title = "original density")
    plot_name = "theoretical_size_vs_in_block_density"
plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")


sparsities = [0.001,0.01,0.03, 0.05,];


plt.figure();
for s in sparsities:
    
    l = 100000;
    a = [0 for i in range(l)];
    
    nz = int(l*s);
    
    for i in range(nz):
        a[i] = 1;
        
    rn.shuffle(a);
    
    x =  [1,5,10,20,30,50,64,80,120, 256];
    y = [];
    y_theory = []
    for block_size in x:
        block_n = int(l/block_size);
    
        print("size", block_size, "n", block_n);
        
        
        full_count = 0;
        for bi in range(block_n):
            for i in range(block_size):
                if a[bi*block_size + i] != 0:
                    full_count += 1;
                    break;
        y.append(full_count/block_n);
        
        yt = ((block_n - 1)/block_n)**(l*s)
        yt = 1 - yt;
        y_theory.append(yt);
        
    
    plt.plot(x,y_theory, label = s)     
    plt.ylabel("percentage of nonzero blocks");
    plt.xlabel("size of blocks");
    
    plt.legend(title = "original density")
    plot_name = "theoretical_size_vs_nonzero_blocks_density"
plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")

