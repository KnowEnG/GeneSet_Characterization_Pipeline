* * * 
## How to verify this pipeline installation on your computer
Use verification testing to assure that the runtime environment and the current version produce the expected output using this repository's data.
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

### 3. Change directory to GeneSet_Characterization_Pipeline/test

```
cd GeneSet_Characterization_Pipeline/test
```

### 4. Start the verification test from the command line

```
make verification_tests
```

### 5. The output files will be compared with the GeneSet_Characterization_Pipeline/data/verification/.../ data
* Each Benchmark will report PASS or FAIL and list the names of files producing differences (if any).
* Note that the files generated will be erased after each Benchmark test.

