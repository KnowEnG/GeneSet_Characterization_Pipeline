#.PHONY: all_unit_tests

all_unit_tests:
	PYTHONPATH='../../src/' python3 -m unittest discover --pattern=test_*.py

clean_outputs:
	rm ../../data/verification/cosine_matrix*
	rm ../../data/verification/DRaWR* 
	rm ../../data/verification/fisher* 
	rm ../../data/verification/net_path* 

