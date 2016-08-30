MKDIR_P = mkdir -p
RMR = rm -rf
GZIP = gzip
GZIP_D = gzip -df
INPUT_DATA_DIR = ./input_data
RUN_DIR = ./run_dir
RESULTS_DIR = $(RUN_DIR)/results
BUILD = ./build
#BUILD = ./test/benchmarks
SCRIPT = ./src/geneset_characterization.py

.PHONY: preparation run_fisher run_drawer decompress_input_data compress_input_data create_run_dir copy_run_files clean_dir_recursively final_clean 

preparation: decompress_input_data create_run_dir copy_run_files

run_fisher:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file fisher_run_file.yml 

run_drawer:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file DRaWR_run_file.yml

decompress_input_data:
	$(GZIP_D) $(INPUT_DATA_DIR)/*
 
compress_input_data:
	$(GZIP) $(INPUT_DATA_DIR)/*

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(BUILD)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

final_clean: clean_dir_recursively compress_input_data
	 
