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
* Create run_paramerters file (yml format) 

 file_name.yml

* Make sure the directories of the input data in `fisher_run_file.yml` and `DRaWR_run_file.yml` are correct
 
* Run GeneSet_Characterization_Pipeline

 ```
 export PYTHONPATH='../GeneSet_Characterization_Pipeline/src':$PYTHONPATH    
 python3 ../GeneSet_Characterization_Pipeline/src/geneset_characterization.py -run_directory ./ -run_file file_name.yml
 ```
  
* Output files are saved in results directory</br>
1.`DRaWR_result` output file saves sorted properties based on the difference between updated user gene vector and baseline.</br>
2.`fisher_result` output file has seven columns and it is sorted in ascending order based on `pval`.

 |**user gene**|**property**|**count**|**user count**|**gene count**|**overlap**|**pval**|
 | ------------| -----------| --------| -------------| -------------| ---------|------|
 |   string    |   string   |    int    |    int    |   int     |   int    |   float    |



### Net_One
