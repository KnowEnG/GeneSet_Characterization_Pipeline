# GeneSet_Characterization_Pipeline
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection

## Steps to run pipelines
1. Configure your environment to have the following packages
2. Run makefile targets
  * Prepare input data and running directories. 
  ```
        make preparation
  ```
  * Run the pipeline you desire
  ```
        make run_fisher
        make run_drawer
  ```
  * Clean the running environment after you finish your tests and analysis
  ```
        make clean_dir_recursively
  ```
 

##DRaWR
##Fisher
##Net_One
