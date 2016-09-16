### How To Generate _html_ Documentation After Cloning the Repository

* Install sphinx, at the command prompt type:
    * >> pip install sphinx
    
* Change to the Samples_Clustering_Pipeline/docs directory in your cloned code.
* Generate the intermediate (.rst) files at the command prompt:
    * >> sphinx-apidoc -f -o source/ ../src/
    
* Use the ./docs/Makefile to generate the _html_ code:
    * >> make html

* In your browser bookmark the ./docs/build/html/index.html documentation homepage.
