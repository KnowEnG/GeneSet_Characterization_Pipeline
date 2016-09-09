# Gene Set Characterization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Gene Set Characterization Pipeline.</br>
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection.

One can select one of three gene set characterization options:

| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Fisher exact test                                | Fisher                               | fisher         |
| Discriminative Random Walks with Restart         | DRaWR                                | DRaWR          |
| Net One                                          | Net One                              | net_one        |

* * * 
## How to run this pipeline with Our data
* * * 
###1. Get Access to KnowEnG-Research Repo:
Email omarsobh@illinois.edu infrastructure team (IST) lead to get:
* __Access__ to KnowEnG-Research github repo

###2. Configurate your Environment in *Two Ways* and Select the One You Prefer
  * Pull the docker image from Dockerhub with the pre-configured environment (__Recommended__)
    1. Install Docker engine in your machine based on your operation system: https://docs.docker.com/engine/installation/
    2. Run the following command to pull our docker image
    
       ```
       docker pull knowengdev/geneset_characterization:09_08_2016
       ```
       
    3. Or install the following package one by one in your machine (Ubuntu or Linux)
  
       ```
       apt-get install -y python3-pip
       pip3 install numpy==1.11.1
       pip3 install pandas==0.18.1 
       pip3 install scipy==0.18.0
       pip3 install scikit-learn==0.17.1
       pip3 install matplotlib==1.4.2
       pip3 install pyyaml
       pip3 install knpackage
       ```
       
        __Note__: *Other dependencies may apply due to different operation system.*

###3. Clone the Samples_Clustering_Pipeline Repo after You are Grant with Access
  ```
   git clone https://github.com/KnowEnG-Research/GeneSet_Characterization_Pipeline.git
  ```

###4. Change directory to the GeneSet_Characterization_Pipeline
  ```
   cd GeneSet_Characterization_Pipeline
  ```

###5. Use the following "make" command to create a local directory "run_dir" below test directory and place all the parameters files in it
  * Change to test directory, prepare input data and running directories
  ```
   cd test
  ```
 
  * Run the pipeline you desire locally
    1. Prepare input data and running directories
      ```
      make env_setup    
      ```
      
    2. Run fisher pipeline</br>
      ```
      make run_fisher
      ```
      
    3. Run DRaWR pipeline</br>
      ```
      make run_drawr
      ```
      
  * Run the pipeline you desire in docker
    1.  Prepare input data and running directories
      ```
      make -f Makefile.docker env_setup
      ```
    
    2. Run fisher pipeline
      ```
      make -f Makefile.docker docker_run_fisher
      ```
   
    3. Run DRaWR pipeline
      ```
      make -f Makefile.docker docker_run_drawr
      ```
 
  * Clean the running environment and revert the local repository to original state after you finish your tests and analysis
    ```
     make final_clean 
    ```

## How to run it with your data 
###6. Setup your run environment
* In the **parent** directory of GeneSet_Characterization_Pipeline directory, create a run directory 

  ```
  mkdir run_directory_name
  ```

* Create results directory to save output files under run directory

  ```
  cd run_directory_name
  mkdir results_directory_name
  ```
* Create run_paramerters file with example below (yml format)
  ```
  cp ../GeneSet_Characterization_Pipeline/data/run_files/fisher_run_file.yml fisher_run_file.yml
  cp ../GeneSet_Characterization_Pipeline/data/run_files/DRaWR_run_file.yml DRaWR_run_file.yml
  ```

 | **Key** | **Value** | **Comments** |
 | ------- | --------- | ------------ |
 | method  | DRaWR or fisher   | Choose DRaWR or fisher as the gene set characterization method |
 | pg_network_name_full_path | directory+pg_network_name |Path and file name of the 4 col property file |
 | gg_network_name_full_path | directory+gg_network_name |Path and file name of the 4 col network file(only needed in DRaWR) |
 | spreadsheet_name_full_path | directory+spreadsheet_name|  Path and file name of user supplied gene sets |
 | results_directory | ./run_dir/results_dir | Directory to save the output files |
 | rwr_max_iterations | 500| Maximum number of iterations without convergence in random walk with restart(only needed in DRaWR) |
 | rwr_convergence_tolerence | 0.0001 | Frobenius norm tolerence of spreadsheet vector in random walk(only needed in DRaWR)|
 | rwr_restart_probability | 0.5 | alpha in `Vn+1 = alpha * N * Vn + (1-alpha) * Vo`(only needed in DRaWR) |
pg_network_name = kegg_pathway_property_gene</br>
gg_network_name = STRING_experimental_gene_gene</br>
spreadsheet_name = ProGENI_rwr20_STExp_GDSC_500.rname.gxc</br>
run_dir = run_directory_name</br>
results_dir = results_directory_name
* Make sure the directories of the input data in `fisher_run_file.yml` and `DRaWR_run_file.yml` are correct
 
* Run GeneSet_Characterization_Pipeline under run_directory_name directory

 ```
 export PYTHONPATH='../GeneSet_Characterization_Pipeline/src':$PYTHONPATH    
 python3 ../GeneSet_Characterization_Pipeline/src/geneset_characterization.py -run_directory ./ -run_file file_name.yml
 ```
  
* Output files are saved in results directory</br>
1.`DRaWR_result` output file saves sorted properties based on the difference between updated user gene vector and baseline.</br>
2.`fisher_result` output file has seven columns and it is sorted in ascending order based on `pval`.

| **user gene** | **property** | **count** | **user count** | **gene count** | **overlap** | **pval** |
|:-------------:|:------------:|:---------:|:--------------:|:--------------:|:-----------:|:--------:|
|   string      |   string     |    int    |    int         |   int          |   int       |   float  |



### Net_One
=======
# GeneSet_Characterization_Pipeline
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection

* User submits significance values (p-values) of all genes.
* User also submits one or more annotations of genes.
* System learns annotations linked to significant genes.
* Probabilistic graphical model.
* Specific example: 

  * p-values are expression-phenotype correlations. 
  * Annotations are eqtls under a tf’s binding sites. 
  * Output: TFs linked to phenotype.
