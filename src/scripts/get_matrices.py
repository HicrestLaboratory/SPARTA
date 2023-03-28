
import ssgetpy
import numpy as np

Ns = (20000,100000)
densities = (0.0001,0.1)
experiment_name = f"suitsparse_N{Ns[0]}_{Ns[1]}_dN{densities[0]}_{densities[1]}"

precision = 1000
matrices = []
for n in range(Ns[0],Ns[1],precision):
    rowbounds = (n,n+precision)
    colbounds = (n,n+precision)
    nzbounds = (int(n*n*densities[0]), int(n*n*densities[1]))
    matrices += ssgetpy.search(rowbounds = rowbounds, colbounds = colbounds, nzbounds = nzbounds)


print(f"found len(matrices) matrices with {Ns[0]} < N,M < {Ns[1]} and {densities[0]} < density < {densities[1]}")
for matrix in matrices:
    matrix.download(extract=True, destpath = f"data/{experiment_name}")