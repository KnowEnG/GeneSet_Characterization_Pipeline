# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT = 		../src/geneset_characterization.py
RUN_DIR = 		./run_dir
RUN_FILES = ../data/run_files
RESULTS_DIR = $(RUN_DIR)/results
INPUT_DATA_DIR = ../data
UNIT_TEST_TEMP_DIR = ./unit/tmp
DATA_DIR = ../data/networks ../data/spreadsheets
INPUT_DATA = $(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*))
INPUT_DATA_SUFFIX = $(suffix $(INPUT_DATA))

RMR = 			rm -rf
GZIP = 			gzip
GZIP_D = 		gzip -df
MKDIR_P = 		mkdir -p 


# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files decompress_input_data

decompress_input_data:
	@echo Decompress input data...
ifeq ($(filter-out .gz, $(INPUT_DATA_SUFFIX)),)
	$(GZIP_D) $(INPUT_DATA)
else
	@echo No gzip file needs to be decompressed.
endif	
 
compress_input_data:
	@echo Compress input data...
ifeq ($(filter .gz, $(INPUT_DATA_SUFFIX)),)
	$(GZIP) $(INPUT_DATA)
else
	@echo Target file in already in gzip format. Aborting compressing input data...
endif

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
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_fisher.yml 

run_drawr:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_2_DRaWR.yml

run_netpath:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_3_net_path.yml

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_fisher_small:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_fisher.yml 

run_drawr_small:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_2_DRaWR.yml

run_netpath_small:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_3_net_path.yml

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
