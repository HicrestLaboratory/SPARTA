Using the utils to analyze the clubs experiments


*Reordering results

#collect reorderings, runs analysis on the reordered and blocked matrices
./collect_all_reord_results 

---Usage: $0 [-s script] [-m matrix_dir] [-r reord_result_dir] [-o output_dir] [-b block_size]




*Multiplication results

#copy and sort the mult files. Organizes them in the root_dir: SPARTA/results/results_2024/mult_data
./rearrange_routine_mult_results
echo "Usage: $0 [-r root_dir] [-u routine] [-x clean folders]"
./utils/rearrange_routine_mult_results.sh -x -r ../../../Downloads/outputs_2024-08-01/SbatchMan/outputs/marzola/ -u spmmcsr

#create csv files
./collect_routine_mult_results 
echo "Usage: $0 [-r root_dir] [-u routine] [-m method]"
./utils/collect_routine_mult_results.sh -r results/results_2024/mult_data/ -u spmmcsr