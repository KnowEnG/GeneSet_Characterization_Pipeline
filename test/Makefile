# Commands that will be used in targets
MKDIR_P = mkdir -p
RMR = rm -rf
GZIP = gzip
GZIP_D = gzip -df
DATA_DIR = ../data/networks ../data/spreadsheets
INPUT_DATA = $(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*))
INPUT_DATA_SUFFIX = $(suffix $(INPUT_DATA)) 
# Directories that will be used in targets
RUN_DIR = ./run_dir
RUN_FILES_DIR = ../data/run_files
RUN_FILE_LOCAL = $(RUN_FILES_DIR)/local
RUN_FILE_DOCKER = $(RUN_FILES_DIR)/docker
RESULTS_DIR = $(RUN_DIR)/results
SCRIPT = ../src/geneset_characterization.py

.PHONY: preparation run_fisher run_drawer decompress_input_data compress_input_data create_run_dir copy_run_files_local clean_dir_recursively final_clean 

# Setup the running environment
env_setup: create_run_dir copy_run_files_local decompress_input_data

# Clean the pipeline after running the pipelines you desire
final_clean: clean_dir_recursively compress_input_data

# Run the fisher pipeline 
run_fisher:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file fisher_run_file.yml 

# Run the DRaWR pipeline
run_drawr:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file DRaWR_run_file.yml

# Decompress input data with .gz suffix
decompress_input_data:
	@echo Decompress input data...
ifeq ($(filter-out .gz, $(INPUT_DATA_SUFFIX)),)
	$(GZIP_D) $(INPUT_DATA)
else
	@echo No gzip file needs to be decompressed.
endif	

# Compress input data to be gzip file
compress_input_data:
	@echo Compress input data...
ifeq ($(filter-out .gxc, $(INPUT_DATA_SUFFIX)),)
	$(GZIP) $(INPUT_DATA)
else
	@echo Illegal file format detected. Aborting compressing input data...
endif	

# Greate run directory
create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

# Copy run_files
copy_run_files_local:
	cp $(RUN_FILE_LOCAL)/*.yml $(RUN_DIR)

copy_run_files_docker:
	cp $(RUN_FILE_DOCKER)/*.yml $(RUN_DIR)

# Remove run directory
clean_dir_recursively:
	$(RMR) $(RUN_DIR)

