# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT = 		../src/geneset_characterization.py
RUN_DIR = 		./run_dir
RUN_FILES = ../data/run_files
RESULTS_DIR = $(RUN_DIR)/results
INPUT_DATA_DIR = ../data
UNIT_TEST_TEMP_DIR = ./unit/tmp

RMR = 			rm -rf
GZIP = 			gzip
GZIP_D = 		gzip -df
MKDIR_P = 		mkdir -p 

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files decompress_input_data

decompress_input_data:
	$(GZIP_D) $(INPUT_DATA_DIR)/networks/* $(INPUT_DATA_DIR)/spreadsheets/*
 
compress_input_data:
	$(GZIP) $(INPUT_DATA_DIR)/networks/* $(INPUT_DATA_DIR)/spreadsheets/*

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)
	
# ----------------------------------------------------------------
# Large Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_fisher:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file fisher.yml 

run_drawr:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file DRaWR.yml

run_netpath:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file net_path.yml

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_fisher_small:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_fisher.yml 

run_drawr_small:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_DRaWR.yml

run_netpath_small:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_net_path.yml

# ----------------------------------------------------------------
# - UNIT TESTS RUN SECTION                                       -
# ----------------------------------------------------------------
unit_test_setup:
	$(MKDIR_P) $(UNIT_TEST_TEMP_DIR)

unit_tests:
	cd unit; make all_unit_tests
# ----------------------------------------------------------------
# - CLEAN UP SECTION                                             -
# ----------------------------------------------------------------	
final_clean: clean_dir_recursively compress_input_data

unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)
