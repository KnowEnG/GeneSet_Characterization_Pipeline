# GeneSet_Characterization_Pipeline
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection

## Steps to run pipelines
###1. Setup github access:
__Access__ KnowEnG-Research github repo

###2. Get a copy of the __GeneSet_Characterization_Pipeline__ code, data
__Run__ the following command to get __GeneSet_Characterization_Pipeline__ repo
```
 git clone https://github.com/KnowEnG-Research/GeneSet_Characterization_Pipeline.git
```
    
###3. Configure your environment to have the following packages
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
   
###4. start in the pipeline repo home directory

```
cd GeneSet_Characterization_Pipeline
```

###5. Prepare 'User Provided Gene Set Spreadsheet' and KnowEnG gene set collections
1. Default user dataset and KnowEnG gene set collections include the following files in __input_data__:   

    | File Name                       | Size    | Col1  | Col2 | Col3 | Col4      | Description                            |
    | --------------------------------| ------- | -------- | ------- | ------- | ------------ | -------------------------------------- |
    | STRING_experimental_gene_gene.gz| 25448   | Gene     | Gene    | Weight | Network Name  | Significant protein interaction dataset|
    | kegg_pathway_property_gene.gz   | 3148460 | Property | Gene    | Weight | Property Name | Pathway propery dataset                |
    
    | File Name           | Format | Description                                    |
    | ------------------- | ------ | ---------------------------------------------- |
    | ProGENI_rwr20_STExp_GDSC_500.rname.gxc.gz| Float  | user spread sheet   |
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
    
  * Clean the running environment and revert the local repository to original state after you finish your tests and analysis
    ```
        make final_clean 
    ```
 

##DRaWR
##Fisher
##Net_One
