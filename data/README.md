#### Example input files are as follows.
- example_train_metadata.txt: the clinical metadata of the training dataset, including experimental group, cohort, and other metadata indices.  
- example_train_abundance.txt: the microbial abundance profile of the training dataset.  
- example_test_metadata.txt: the clinical metadata of the external test dataset.  
- example_test_abundance.txt: the microbial abundance profile of the external test dataset.  
- example_other_metadata.txt: the clinical metadata of samples with other diseases used to assess the specificity of identified markers.  
- example_other_abundance.txt: the microbial abundance profile of samples with other diseases.
- covariates.txt: the example input covariates file for Step 4.
- convert.tsv: the example converted microbial count matrices for Step 11.

#### Example results are as follows.
- differential_plot.pdf: the volcano plot of differential signatures (output file of Step 4).  
- cross_validation_auc.pdf: the visualization of the cross-validation AUC of the best-performing model constructed in Step 7.
- validation_AUC.pdf: the visualization of performance metrics of markers in internal validations (output file of Step 8).
- marker_importance.pdf: the visualization of permutation feature importance of candidate markers in the best-performing model (output file of Step 10).
