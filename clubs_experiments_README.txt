Using the utils to analyze the clubs experiments


*Reordering results

#collect reorderings, runs analysis on the reordered and blocked matrices
./collect_all_reord_results 

---Usage: $0 [-s script] [-m matrix_dir] [-r reord_result_dir] [-o output_dir] [-b block_size]




*Multiplication results

#copy and sort the mult files. Organizes them in the root_dir: SPARTA/results/results_2024/mult_data
./rearrange_mult_results

./collect_all_mult_results
echo "Usage: $0 [-r root_dir] [-u routine] [-m method]"
