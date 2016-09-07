# Commands that will be used in targets
MKDIR_P = mkdir -p
RMR = rm -rf
GZIP = gzip
GZIP_D = gzip -df
DATA_DIR = ../data/networks ../data/spreadsheets
OBJS = $(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*))
#OBJS = $(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*))
# Directories that will be used in targets
RUN_DIR = ./run_dir
RUN_FILES = ../data/run_files
RESULTS_DIR = $(RUN_DIR)/results
SCRIPT = ../src/geneset_characterization.py

.PHONY: preparation run_fisher run_drawer decompress_input_data compress_input_data create_run_dir copy_run_files clean_dir_recursively final_clean 

setup_env: decompress_input_data create_run_dir copy_run_files

run_fisher:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file fisher_run_file.yml 

run_drawr:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file DRaWR_run_file.yml

decompress_input_data: $(OBJS)
	if ls -s $(OBJS)  ; \
	then \
		$(GZIP_D) $(OBJS); \
	fi;

compress_input_data:
	ifneq ( $(suffix $(OBJS)), .gz )
		$(GZIP) $(OBJS)
	endif


create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

final_clean: clean_dir_recursively compress_input_data
	 
