#.PHONY: all_unit_tests

all_unit_tests:
	PYTHONPATH='../../src/' python3 -m unittest discover --pattern=test_*.py

clean_outputs:
	rm unit_test_run_dir/results/cosine_matrix.df
	rm unit_test_run_dir/results/DRaWR_droplist* unit_test_run_dir/results/DRaWR_result*
	rm unit_test_run_dir/results/fisher_droplist* unit_test_run_dir/results/fisher_result*
	rm unit_test_run_dir/results/net_path_droplist* unit_test_run_dir/results/net_path_result*

