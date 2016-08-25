# GeneSet_Characterization_Pipeline
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection

## Steps to run pipelines
1. Configure your environment to have the following packages
  ```
    System: ubuntu:14.04
    apt-get install -y python3-pip
    apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
    pip3 install -I numpy==1.11.1
    pip3 install -I pandas==0.18.1 
    pip3 install -I scipy==0.18.0
    pip3 install -I scikit-learn==0.17.1
    pip3 install -I vcversioner==2.16.0.0
    pip3 install -I knpackage==0.1.2
    apt-get install -y libfreetype6-dev libxft-dev 
    pip3 install -I matplotlib==1.4.2
    pip3 install pyyaml
   ```
   
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
