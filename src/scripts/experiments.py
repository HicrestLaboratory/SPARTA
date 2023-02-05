# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import csr as csr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

def evaluate_blocking(cmat, grouping, block_size):
    groups = np.unique(grouping)
    block_count = 0;
    total_block_height = 0;
    total_fill_in = 0;
    for group in groups:
        pattern = []
        group_size = 0

        #build the pattern;
        for row, g in enumerate(grouping):
            if g == group:
                pattern = sims.merge_patterns(pattern, cmat.pos[row])
                group_size += 1
        
                
        #count nz_blocks; add their size
        block_pattern = sims.get_block_pattern(pattern, block_size)
        block_count += len(block_pattern);
        total_block_height += len(block_pattern)*group_size
        
        #count nonzeros
        for row, g in enumerate(grouping):
            if g == group:
                #now count nonzeros
                blocked_row = sims.get_block_pattern(cmat.pos[row],block_size)
                total_fill_in += block_size*len(np.setdiff1d(block_pattern, blocked_row, True));
                total_fill_in += len(blocked_row)*block_size - len(cmat.pos[row])
                
                """
                if blocked_row:
                    if blocked_row[-1] == int(cmat.M/block_size):
                        total_fill_in -= block_size - (cmat.M)%block_size
                        
                if block_pattern:
                    if block_pattern[-1] == int(cmat.M/block_size):
                        total_fill_in -= block_size - (cmat.M)%block_size
                """
        
    return block_count, total_block_height, total_fill_in



def run_experiment(cmat, taus, block_size, name, blocking, dist_func):
    densities = []
    avg_sizes = []    
        
    
    blocked_dist_func = lambda p1,p2,s1,s2 : dist_func(p1,p2,s1,s2,block_size)
        
    for tau in taus:
        #print(f"evaluating for tau {tau}")
        grouping = blocking(cmat, tau, blocked_dist_func)

        block_count, total_block_height, total_fill_in = evaluate_blocking(cmat, grouping, block_size)
        
        density = cmat.tot_nz()/(cmat.tot_nz() + total_fill_in)
        densities.append(density)
        avg_sizes.append(total_block_height/block_count)
    return densities,avg_sizes


IC = {
"name" : "IC-simple",
"blocking" : blk.IterativeClusteringSimple,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.JaccardGeneral(p1,p2),
}

IC_quotient = {
"name" : "IC quotient",
"blocking" : blk.IterativeClusteringSimple,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.JaccardGeneral(p1,p2,block_size = block_size),
}


ICP = {
"name" : "IC",
"blocking" : blk.IterativeClusteringPattern,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.JaccardGeneral(p1,p2,block_size = block_size)
}

ICP_special = {
"name" : "IC pattern special",
"blocking" : blk.IterativeClusteringPattern,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.JaccardGeneral(p1,p2,block_size = block_size, s1 = s1)
}

HC = {
"name" : "HC-simple",
"blocking" : blk.HierarchicalClustering,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.HammingGeneral(p1,p2,block_size = block_size)
}

HCH = {
"name" : "HC",
"blocking" : blk.HierarchicalClustering,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.HammingGeneral(p1,p2,block_size = block_size, s1 = s1, s2 = s2)
}


HCJ = {
"name" : "HC - Jaccard",
"blocking" : blk.HierarchicalClustering,
"dist_func" : lambda p1,p2,s1,s2,block_size : sims.JaccardGeneral(p1,p2,block_size = block_size, s1 = s1, s2 = s2)
}



def BCSR_blocking(cmat, block_height):
    ga = []
    i = 0
    while (len(ga) < cmat.N):
        ga += [i,]*block_height
        i += 1
    return ga[0:cmat.N]

        
    

def make_image(m,n,mat_density, block_size = 16, tau_amount = 50):
    
    colors = itertools.cycle(sns.color_palette('deep'))
    
    cmat = csr.make_random_CSR(m,n,mat_density);
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    for exp in [HCH, ICP, IC]:
        taus = np.linspace(0, 1, tau_amount);
        if "HC" in exp["name"]:
            taus = np.linspace(0,cmat.M,50)
        densities, sizes = run_experiment(cmat, taus, block_size, **exp)
        ax.plot(densities, sizes, label = exp["name"], color = next(colors))
        print(f"processed {exp['name']}")
        
    #add BCSR experiment
    densities = []
    avg_sizes = [] 
    
    for block_height in np.linspace(1, cmat.N, cmat.N):
        block_height = int(block_height)    
        grouping = BCSR_blocking(cmat, block_height)
        block_count, total_block_height, total_fill_in = evaluate_blocking(cmat, grouping, block_size)
        density = cmat.tot_nz()/(cmat.tot_nz() + total_fill_in)
        densities.append(density)
        avg_sizes.append(total_block_height/block_count)
        
        
    ax.plot(densities, avg_sizes, label = "FIXED", color = next(colors))
    
    
    
    #plt.figure(figsize=(8,4), tight_layout=True)
    ax.set_yscale("log")
    ax.set_xlabel("density inside nonzero blocks")
    ax.set_ylabel("avg. size of nonzero blocks")
    ax.legend()
    plt.savefig(f'reordering_curves_{m}x{n}_{mat_density}_{block_size}.png', dpi = 300)    






def make_image_from_edgelist(filename, block_size = 16, tau_amount = 50):
    
    colors = itertools.cycle(sns.color_palette('deep'))
    
    cmat = csr.CSR()
    cmat.fill_from_edgelist(filename, delimiter = " ")
    
    taus = np.linspace(0, 1, tau_amount);
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    for exp in [IC, ICP]:
        densities, sizes = run_experiment(cmat, taus, block_size, **exp)
        ax.plot(densities, sizes, label = exp["name"], color = next(colors))
        print(f"processed {exp['name']}")
        
    #add BCSR experiment
    densities = []
    avg_sizes = [] 
    
    for block_height in np.linspace(1, cmat.N, cmat.N):
        block_height = int(block_height)    
        grouping = BCSR_blocking(cmat, block_height)
        block_count, total_block_height, total_fill_in = evaluate_blocking(cmat, grouping, block_size)
        density = cmat.tot_nz()/(cmat.tot_nz() + total_fill_in)
        densities.append(density)
        avg_sizes.append(total_block_height/block_count)
        
        
    ax.plot(densities, avg_sizes, label = "FIXED", color = next(colors))
    
    
    
    #plt.figure(figsize=(8,4), tight_layout=True)
    ax.set_yscale("log")
    ax.set_xlabel("density inside nonzero blocks")
    ax.set_ylabel("avg. size of nonzero blocks")
    ax.legend()
    plt.savefig(f'reordering_curves_{filename.split(".")[0]}_{block_size}.png', dpi = 300)    




"""
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
plt.rc('legend', fontsize=13)    # legend fontsize
plt.rc('font', size=13)          # controls default text sizes
"""
#make_image(128,128,0.1,8)
#make_image(256,256,0.05,8)
make_image(256,256,0.05,8)
make_image(512,512,0.01,32)

#make_image(256,256,0.05,16)
#make_image(256,256,0.05,8)

#make_image(256,256,0.1,16)
#make_image(256,256,0.1,8)


#make_image(512,512,0.01,32)
#make_image(512,512,0.05,64)
#make_image(1024,1024,0.05,64)

make_image_from_edgelist("ia-wikiquote-user-edits-nodup.el", block_size = 1)

