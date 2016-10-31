# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT = 		../src/geneset_characterization.py
RUN_DIR = 		./run_dir
DATA_DIR = 		../data/networks ../data/spreadsheets
INPUT_DATA = 		$(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*))
RESULTS_DIR = 		$(RUN_DIR)/results
RUN_FILES_DIR = 	../data/run_files
RUN_FILE_LOCAL = 	$(RUN_FILES_DIR)
RUN_FILE_DOCKER = 	$(RUN_FILES_DIR)/docker
INPUT_DATA_SUFFIX = 	$(suffix $(INPUT_DATA)) 

RMR = 			rm -rf
GZIP = 			gzip
GZIP_D = 		gzip -df
MKDIR_P = 		mkdir -p

.PHONY: preparation run_fisher run_drawer decompress_input_data compress_input_data create_run_dir copy_run_files_local clean_dir_recursively final_clean 

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files_local decompress_input_data

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
	@echo Illegal file format detected. Aborting compressing input data...
endif

copy_run_files_local:
	cp $(RUN_FILE_LOCAL)/*.yml $(RUN_DIR)

copy_run_files_docker:
	cp $(RUN_FILE_DOCKER)/*.yml $(RUN_DIR)

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 
	
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
# - CLEAN UP SECTION                                             -
# ----------------------------------------------------------------	
final_clean: clean_dir_recursively compress_input_data

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

