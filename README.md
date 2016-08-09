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
