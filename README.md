# KnowEnG's Gene Set Characterization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Gene Set Characterization Pipeline.

This pipeline **ranks** a user supplied gene set **against** a KnowEnG's gene sets collection.

There are three gene set characterization methods that one can choose from:

| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Fisher exact test                                | Fisher                               | fisher         |
| Discriminative Random Walks with Restart         | DRaWR                                | DRaWR          |
| Net One                                          | Net One                              | net_one        |

* * * 
## How to run this pipeline with Our data
* * * 
###1. Get Access to KnowEnG-Research Repo:
Email omarsobh@illinois.edu infrastructure team (IST) lead to:

* __Access__ KnowEnG-Research github repo

###2. Clone the GeneSet_Characterization_Pipeline Repo
```
 git clone https://github.com/KnowEnG-Research/GeneSet_Characterization_Pipeline.git
```
 
###3. Install the following (Ubuntu or Linux)
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

###4. Change directory to GeneSet_Characterization_Pipeline

```
cd GeneSet_Characterization_Pipeline
```

###5. Change directory to test

```
cd test
```
 
###6. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

###7. Select and run a gene set characterization option:
 
 * Run fisher pipeline</br>
  ```
  make run_fisher
  ```
 
 * Run DRaWR pipeline</br>
  ```
  make run_drawr
  ```


* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-4 above then do the following:***__

### * Create your run directory

 ```
 mkdir run_directory
 ```

### * Change directory to the run_directory

 ```
 cd run_directory
 ```

### * Create your results directory

 ```
 mkdir results_directory
 ```
 
### * Create run_paramters file (YAML Format)
 ``` 
 Look for examples of run_parameters in the GeneSet_Characterization_Pipeline/data/run_files template_run_parameters.yml
 ```

### * Run the GeneSet Characterization Pipeline:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='./src':$PYTHONPATH    
   ```
   
  * Run
   ```
  python3 ./src/geneset_characterization.py -run_directory ./ -run_file template_run_parameters.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                   | **Value** | **Comments** |
| ------------------------- | --------- | ------------ |
| method                    | DRaWR or fisher   | Choose DRaWR or fisher as the gene set characterization method |
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
results_dir = results_directory_name</br>

* * * 
## Description of Output files saved in results directory
* * * 

* `DRaWR_result` output file saves sorted properties based on the difference between updated user gene vector and baseline.</br>

 | **user gene set name1** |**user gene set name2**|**...**|**user gene set name n**|**base**|
 | :--------------------: |:--------------------:|---|:--------------------:|:--------------------:|
 | property name (string)</br> (most significant) |property name (string)</br> (most significant)|...|property name (string)</br> (most significant)|property name (string)</br> (most significant)|
 | ... |...|...|...|...|
 | property name (string)</br> (least significant) |property name (string)</br> (least significant)|...|property name (string)</br> (least significant)|property name (string)</br> (least significant)|
* `fisher_result` output file has seven columns and it is sorted in ascending order based on `pval`.

 | **user gene** | **property** | **count** | **user count** | **gene count** | **overlap** | **pval** |
 |:-------------:|:------------:|:---------:|:--------------:|:--------------:|:-----------:|:--------:|
 |   string      |   string     |    int    |    int         |   int          |   int       |   float  |
