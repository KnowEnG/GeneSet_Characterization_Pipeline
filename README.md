# KnowEnG's Gene Set Characterization Pipeline

This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Gene Set Characterization Pipeline.

This pipeline **ranks** a user supplied gene set **against** a KnowEnG's gene sets collection.

There are three gene set characterization methods that one can choose from:

| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Fisher exact test                                | Fisher                               | fisher         |
| Discriminative Random Walks with Restart         | DRaWR                                | DRaWR          |
| Net Path                                         | Net Path                             | net_path       |

* * * 
## How to run this pipeline with Our data
* * * 
 
### 1. Clone the GeneSet_Characterization_Pipeline Repo
```
 git clone https://github.com/KnowEnG/GeneSet_Characterization_Pipeline.git
```
### 2. Install the following (Ubuntu or Linux)
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

### 3. Change directory to GeneSet_Characterization_Pipeline

```
cd GeneSet_Characterization_Pipeline
```

### 4. Change directory to test

```
cd test
```
 
### 5. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

### 6. Select and run a gene set characterization option:
 
 * Run fisher pipeline</br>
  ```
  make run_fisher
  ```
 
 * Run DRaWR pipeline</br>
  ```
  make run_drawr
  ```
  
 * Run Net Path pipeline</br>
  ```
  make run_netpath
  ```


* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-3 above then do the following:***__

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
 Look for examples of run_parameters in the GeneSet_Characterization_Pipeline/data/run_files BENCHMARK_1_fisher.yml
 ```

### * Run the GeneSet Characterization Pipeline:

  * Update PYTHONPATH environment variable
   ``` 
   export PYTHONPATH='./src':$PYTHONPATH    
   ```
   
  * Run
   ```
  python3 ../src/geneset_characterization.py -run_directory ./run_dir -run_file BENCHMARK_1_fisher.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                   | **Value** | **Comments** |
| ------------------------- | --------- | ------------ |
| method                    | DRaWR or fisher or net_path   | Choose DRaWR or fisher or Net Path as the gene set characterization method |
| pg_network_name_full_path | directory+pg_network_name |Path and file name of the 4 col property file |
| gg_network_name_full_path | directory+gg_network_name |Path and file name of the 4 col network file(needed in DRaWR and Net Path) |
| spreadsheet_name_full_path | directory+spreadsheet_name|  Path and file name of user supplied gene sets |
| gene_names_map | directory+gene_names_map| Map ENSEMBL names to user specified gene names |
| results_directory | directory | Directory to save the output files |
| rwr_max_iterations | 500| Maximum number of iterations without convergence in random walk with restart(needed in DRaWR and Net Path) |
| rwr_convergence_tolerence | 0.0001 | Frobenius norm tolerence of spreadsheet vector in random walk(needed in DRaWR and Net Path)|
| rwr_restart_probability | 0.5 | alpha in `V_(n+1) = alpha * N * Vn + (1-alpha) * Vo` (needed in DRaWR and Net Path) |
| k_space| 100| number of the new space dimensions in SVD(only needed in Net Path)|

pg_network_name = kegg_pathway_property_gene.edge</br>
gg_network_name = STRING_experimental_gene_gene.edge</br>
spreadsheet_name = ProGENI_rwr20_STExp_GDSC_500.rname.gxc.tsv</br>
gene_names_map = ProGENI_rwr20_STExp_GDSC_500_MAP.rname.gxc.tsv

* * * 
## Description of Output files saved in results directory
* * * 

* Output files of all three methods save sorted properties for each gene set with name {method}\_ranked\_by\_property\_{timestamp}.df.</br>

 | **user gene set name1** |**user gene set name2**|**...**|**user gene set name n**|
 | :--------------------: |:--------------------:|---|:--------------------:|
 | property </br> (most significant) |property </br> (most significant)|...|property </br> (most significant)|
 | ... |...|...|...|
 | property</br> (least significant) |property </br> (least significant)|...|property</br> (least significant)|
* Fisher method saves one output file with seven columns and it is sorted in descending order based on `pval`. The name of the file is fisher_sorted_by_property_score\_{timestamp}.df. 

 | **user_gene_set** | **property_gene_set** | **pval** | **universe_count** | **user_count** | **property_count** | **overlap_count** |
 |:-------------:|:------------:|:---------:|:--------------:|:--------------:|:-----------:|:--------:|
 |   user gene 1      |   property 1    |    float    |    int         |   int          |   int       |   int  |
 |    ...             |   ...           |    ...      |    ...         |   ...          |  ...        |   ...  | 
 |   user gene n      |   property m    |    float    |    int         |   int          |   int       |   int  |
 
* DRaWR method saves two output files with five columns and they are sorted in descending order based on `difference_score`. The files are DRaWR_sorted_by_gene_score\_{timestamp}.df and DRaWR_sorted_by_property_score\_{timestamp}.df

 | **user_gene_set** | **gene_node_id** | **difference_score** | **query_score** | **baseline_score** |
 |:-------------:|:------------:|:---------:|:--------------:|:--------------:|
 |   user gene 1      |   gene node 1     |    float    |    float         |   float          | 
 |    ...      |   ...     |    ...    |    ...         |   ...          | 
 |   user gene n      |   gene node m     |    float    |    float         |   float          | 
 
 | **user_gene_set** | **property_gene_set** | **difference_score** | **query_score** | **baseline_score** |
 |:-------------:|:------------:|:---------:|:--------------:|:--------------:|
 |   user gene 1      |   property 1     |    float    |    float         |   float          | 
 |    ...      |   ...     |    ...    |    ...         |   ...          | 
 |   user gene n      |   property m     |    float    |    float         |   float          | 
 
 
* Net Path method saves one output file with three columns and it is sorted in descending order based on `cosine_sum`. The name of the file is net_path_sorted_by_property_score\_{timestamp}.df. 

 | **user_gene_set** | **property_gene_set** | **cosine_sum** |
 |:-------------:|:------------:|:---------:| 
 |   user gene 1      |   property 1     |    float    | 
 |...|...|...|
 |user gene n| property m| float|
 
 
 
