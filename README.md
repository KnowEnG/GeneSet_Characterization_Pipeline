<<<<<<< HEAD
# Gene Set Characterization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Gene Set Characterization Pipeline.</br>
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection.

One can select one of three gene set characterization options:

| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Fisher exact test                                | Fisher                               | fisher         |
| Discriminative Random Walks with Restart         | DRaWR                                | DRaWR          |
| Net One                                          | Net One                              | net_one        |

## How to run this pipeline with provided data
###1. Get Access to KnowEnG-Research Repo:
Email omarsobh@illinois.edu infrastructure team (IST) lead to:

1. __Access__ KnowEnG-Research github repo

###2. Clone the Samples_Clustering_Pipeline Repo
```
 git clone https://github.com/KnowEnG-Research/GeneSet_Characterization_Pipeline.git
```
###3. Install the following (Mac OS or Linux)
```
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1 
 pip3 install scipy==0.18.0
 pip3 install scikit-learn==0.17.1
 apt-get install -y libfreetype6-dev libxft-dev 
 pip3 install matplotlib==1.4.2
 pip3 install pyyaml
 pip3 install knpackage
```

###4. Change directory to  the GeneSet_Characterization_Pipeline
```
 cd GeneSet_Characterization_Pipeline
```

###5. Use the following "make" command to create a local directory "run_dir" and place all the parameters files in it
  * Prepare input data and running directories. 
 ```
  make preparation
 ```
 
  * Run the pipeline you desire
 ```
  make -f Makefile.docker docker_run_fisher
  make -f Makefile.docker docker_run_drawer
 ```
 
  * Clean the running environment and revert the local repository to original state after you finish your tests and analysis
 ```
  make final_clean 
 ```

## How to run it with your data 
###6. Setup your run environment
* Create a  run directory

 ```
 mkdir run_directory_name
 ```

* Create results directory to save output files under run directory

 ```
 cd run_directory_name
 mkdir results_directory_name
 ```
* Create run_paramerters file with example below (yml format)
 
 | **Key** | **Value** | **Comments** |
 | ------- | --------- | ------------ |
 | method  | DRaWR or fisher   | Choose DRaWR or fisher as the gene set characterization method |
 | pg_network_name_full_path | ./a/pg_network_name |Path and file name of the 4 col property file |
 | gg_network_name_full_path | ./a/gg_network_name |Path and file name of the 4 col network file(only needed in DRaWR) |
 | spreadsheet_name_full_path | ./a/spreadsheet_name|  Path and file name of user supplied gene sets |
 | results_directory | ./run_dir/results | Directory to save the output files |
 | rwr_max_iterations | 500| Maximum number of iterations without convergence in random walk with restart(only needed in DRaWR) |
 | rwr_convergence_tolerence | 0.0001 | Frobenius norm tolerence of spreadsheet vector in random walk(only needed in DRaWR)|
 | rwr_restart_probability | 0.5 | alpha in `Vn+1 = alpha * N * Vn + (1-alpha) * Vo`(only needed in DRaWR) |
a = input_date</br>
pg_network_name = kegg_pathway_property_gene</br>
gg_network_name = STRING_experimental_gene_gene</br>
spreadsheet_name = ProGENI_rwr20_STExp_GDSC_500.rname.gxc
* Make sure the directories of the input data in `fisher_run_file.yml` and `DRaWR_run_file.yml` are correct
 
* Run GeneSet_Characterization_Pipeline

 ```
 export PYTHONPATH='../GeneSet_Characterization_Pipeline/src':$PYTHONPATH    
 python3 ../GeneSet_Characterization_Pipeline/src/geneset_characterization.py -run_directory ./ -run_file file_name.yml
 ```
  
* Output files are saved in results directory</br>
1.`DRaWR_result` output file saves sorted properties based on the difference between updated user gene vector and baseline.</br>
2.`fisher_result` output file has seven columns and it is sorted in ascending order based on `pval`.

 | **user gene** | **property** | **count** | **user count** | **gene count** | **overlap** | **pval** |
 | ------------- | ------------ | --------- | -------------- | -------------- |  ---------- | -------- |
 |   string    |   string   |    int    |    int    |   int     |   int    |   float    |



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
>>>>>>> upstream/master
