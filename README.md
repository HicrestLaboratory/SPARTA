Compile with 

'''make serial''' to compile without cuda  
or
'''make all''' to compile also the cuda test

run '''./programs/general/TEST_blocking_VBR''' to see an example application; you can see available options in input.h

important options: 
-a algorithm choice (0: iterative, 1: structured, 2: fixed)
-b column_block_size
-m similarity_function (0:Hamming 1:Jaccard 2:HammingOPENMP 3:JaccardOPENMP
