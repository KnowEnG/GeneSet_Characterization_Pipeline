#.PHONY: all_unit_tests

all_unit_tests:
	PYTHONPATH='../../src/' python3 -m unittest discover --pattern=test_*.py

clean_outputs:
	rm DRaWR_result* fisher_result* 
	rm unit_test_run_dir/DRaWR_droplist* unit_test_run_dir/DRaWR_result* unit_test_run_dir/fisher_droplist* unit_test_run_dir/fisher_result*

