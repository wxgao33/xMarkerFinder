# Identification and validation of microbial marker from cross-cohort datasets using xMarkerFinder
xMarkerFinder is a four-stage workflow for microbiome research including cross-cohort differential signature analysis, model construction, model validation, and marker identification.  



## User Tutorial
### Stage 1 Cross-cohort differential signature analysis  
#### 1. Data normalization. 
Convert microbial counts to relative abundance profiles of all datasets involved.  
```
$ Rscript 1_Normalization.R -W /workplace/ -p abundance.txt -o TEST
```  
Users should specify these parameters or enter the default values, subsequent repetitions of which are not listed.   
```
-W the Workplace of this whole protocol  
-p the input microbial count profile  
-o prefix of output files
```  
- Input files:  
abundance.txt: merged microbial count profile of all datasets.  
- Output files:  
relative_abundance.txt: normalized relative abundance profile of input dataset. Relative abundance profiles are used as input files for all subsequent analyses, except for Steps 20-21, which requires raw count file.  
#### 2.	Data filtering. 
Filter microbial signatures with low occurrence rates across cohorts.  
```
$ Rscript 2_Filtering.R -W /workplace/ -m metadata.txt -p relative_abundance.txt -b Cohort -t 2 -o TEST  
```
```
-m the input metadata file  
-p the input microbial relative abundance file (output file of Step 1)  
-b the column name of batch(cohort) in metadata (default: Cohort)  
-t the minimum number of cohorts where features have to occur (default: 2)  
-O prefix of output files  
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
relative_abundance.txt: normalized relative abundance profile of the training dataset.  
- Output files:  
filtered_abundance.txt: filtered relative abundance profile of the training dataset, used as the input file for following steps.  
#### 3.	Confounder analysis.   
PERMANOVA test based on Bray-Curtis dissimilarity is performed to assess the correlation between metadata and microbial profiles and returns the coefficient of determination (R2) value and P value of every metadata index. Whichever index that contributes the most is considered as the major confounder and used later in Step 4. PCoA plot with Bray-Curtis dissimilarity is provided.  
```
$ Rscript 3_Confounder_analysis.R -W /workplace/ -m metadata.txt -p filtered_abundance.txt -g Group -o TEST  
```
```
-m input metadata file  
-p input filtered microbial relative abundance file  
-g the column name of experimental interest(group) in metadata (default: Group)  
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
filtered_abundance.txt: filtered relative abundance profile after preprocessing.  
- Output files:  
metadata_microbiota.txt: the confounding effects caused by clinical information, used to determine the major batch and covariates.  
#### 4.	Differential analysis.   
Based on the major confounder and covariates found in Step 3, cross-cohort differential signature analysis is conducted.
```
$ Rscript 4_Differential_analysis.R -W /workplace/ -m metadata.txt -p filtered_abundance.txt -g Group -b Cohort -c covariates.txt -o TEST
```
```
-g the column name of experimental interest(group) in metadata (default: Group)  
-b the column name of major confounder in metadata (default: Cohort)  
-c input covariates file (tab-delimited format containing all covariates)  
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
filtered_abundance.txt: filtered relative abundance profile after preprocessing.  
covariates.txt: covariates identified in Step 3 (newly generated tab-delimited file where each row is a covariate, example file is provided).  
- Output files:  
differential_significance_single_cohort.txt: the differential significance result in individual cohorts.  
differential_significance.txt: meta-analytic testing results aggregating differential testing results in individual cohorts, used for next-step visualization.  
differential_signature.txt: significantly differential signatures between groups derived from input filtered profiles, used as input files for feature selection.  
#### 5.	Differential signature volcano plot.  
This step provides the visualization of Step 4 with a volcano plot where the identified differential signatures are marked as red (upregulated in case group) and blue (downregulated in case group) dots.
```
$ python 5_Differential_signature_volcano_plot.py -W /workplace/ -i differential_significance.txt -t 0.05 -o TEST
```
```
-i the input differential significance file (output file of Step 4)
-t the threshold of P value for plotting (default: 0.05)
```
- Input files:  
differential_significance.txt: meta-analytic differential testing results.  
- Output files:  
differential_volcano.pdf: the volcano plot of input differential significance file.  
### Stage 2 Model construction
#### 6.	Classifier selection.   
This step provides optional classifier selection for subsequent steps where the performances of every ML algorithm are generally assessed using all differential signatures. The output file contains the cross-validation AUC, specificity, sensitivity, accuracy, precision and F1 score of all classification models built with these various algorithms. Users should specify the selected classifier in all following steps.
```
$ python 6_Classifier_selection.py -W /workplace/ -m metadata.txt -p differential_signature.txt -g Group -e exposure -s 0 -o TEST
```
```
-p input differential signature file (output file of Step 4)
-g the column name of experimental interest(group) in metadata (default: Group)
-e the experiment group(exposure) of interest (in example data: CRC)
-s random seed (default:0)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
differential_signature.txt: significantly differential signatures between groups.  
- Output files:  
classifier_selection.txt: the overall cross-validation performance of all classifiers using differential signatures, used to determine the most suitable classifier.  
#### 7.	Feature effectiveness evaluation.
The first step of Triple-E feature selection procedure evaluates the predictive capability of every feature via constructing individual classification models respectively. Users should specify an ML algorithm here and in every future step as the overall classifier for the whole protocol from the following options: Logisticl1, Logisticl2, DecisionTree, RandomForest, GradientBoost, KNeighbors and SVC. Features with cross-validation AUC above the threshold (default:0.5) are defined as effective features and are returned in the output file.  
```
$ python 7_Feature_effectiveness_evaluation.py -W /workplace/ -m metadata.txt -p differential_signature.txt -g Group -e exposure -b Cohort -c classifier -s 0 -t 0.5 -o TEST
```
```
-p input differential signature file (output file of Step 4)
-b the column name of batch(cohort) in metadata (default: Cohort)
-c selected classifier
-t AUC threshold for defining if a feature is capable of prediction (default:0.5)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
differential_signature.txt: significantly differential signatures between groups.  
- Output files:  
feature_auc.txt: cross-validation AUC values of individual features.  
effective_feature.txt: features derived from differential signatures that are capable of predicting disease states, used as input file of the following step.  
#### 8.	Collinear feature exclusion.   
The second step of feature selection aims to exclude collinear issue caused by highly correlated features based on the result of Step 7 and returns the uncorrelated-effective features.
```
$ python 8_Collinear_feature_exclusion.py -W /workplace/ -p effective_feature.txt -t 0.7 -o TEST
```
```
-p input effective feature file (output file of Step 7)
-t correlation threshold for collinear feature exclusion (default:0.7)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
effective_feature.txt: features with classification capability.  
- Output files:  
feature_correlation.txt: spearman correlation coefficient of every feature pair.  
uncorrelated_effective_feature.txt: features derived from input effective features excluding highly collinear features, used as input file of the following step.  
#### 9.	Recursive feature elimination.   
The last step of feature selection recursively eliminates the weakest feature per loop to sort out the minimal panel of candidate markers.  
```
$ python 9_Recursive_feature_elimination.py -W /workplace/ -m metadata.txt -p uncorrelated_effective_feature.txt -g Group -e exposure -c classifier -s 0 -o TEST
```
```
-p input uncorrelated-effective feature file (output file of Step 8)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
uncorrelated_effective_feature.txt: independent features derived from effective features.  
- Output files:  
candidate_marker.txt: identified optimal panel of candidate markers, used as model input for all subsequent steps.  
#### 10.	Boruta feature selection. 
Besides Triple-E feature selection procedure, we provide an alternative method, feature selection with the Boruta algorithm.  
```
$ Rscript 10_Boruta_feature_selection.R -W /workplace/ -m metadata.txt -p differential_signature.txt -g Group -s 0 -o TEST
```
```
-p input differential signature profile (output file of Step 4) or uncorrelated-effective feature file (output file of Step 8)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
differential_signature.txt: differential signatures used for feature selection (could also be uncorrelated-effective features from Step 8).  
- Output files:  
boruta_feature_imp.txt: confirmed feature importances via Boruta algorithm.  
boruta_selected_feature.txt: selected feature profile, used as input candidate markers for subsequent steps.  
#### 11.	Hyperparameter tuning and cross-validation of the best-performing model.   
Based on the selected classifier and candidate markers, the hyperparameters of the classification model are adjusted via bayesian optimization method based on cross-validation AUC. The output files contain the tuned hyperparameters and the multiple performance metric values of the constructed best-performing model.  
```
$ python 11_Model_construction.py -W /workplace/ -m metadata.txt -p candidate_marker.txt -g Group -e exposure -c classifier -s 0 -o TEST
```
```
-p input candidate marker profile (output file of Step 9 or Step 10)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers (or boruta_selected_feature.txt for all subsequent steps).  
- Output files:  
best_param.txt: the best hyperparameter combination of classification model.  
optimal_cross_validation.txt: the overall cross-validation performance of the best-performing model.  
#### 12.	Cross validation AUC plot. 
This step provides the visualization for the cross-validation AUC of the best-performing model constructed in Step 11.  
```
$ python 12_Cross_validation_AUC_plot.py -W /workplace/ -m metadata.txt -p candidate_marker.txt -g Group -e exposure -c classifier -r best_param.txt -s 0 -o TEST
```
```
-p input candidate feature file (output file of Step 9 or Step 10)
-r input optimal hyperparameter file (output file of Step 11)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers.  
best_param.txt: the best hyperparameter combination of classification model.  
- Output files:  
cross_validation_auc.pdf: the visualization of the cross-validation AUC of the best-performing model.  







## FAQs
### Part I General questions
#### 1. When should I use xMarkerFinder?  
xMarkerFinder is suitable for microbial marker identification from cross-cohort datasets. Our previous studies demonstrated its applicability in identifying global microbial diagnostic markers for adenoma and colorectal cancer. Moreover, xMarkerFinder could also be applied to marker determination in disease prognosis, treatment stratification, metastasis surveillance, adverse reactions anticipation, etc. Any research dedicated to marker identification from multi-population microbial datasets are welcome.
#### 2. How should I setup the required computational environment for xMarkerFinder?  
We provide detailed instructions on software installation for users to run the whole xMarkerFinder workflow locally. However, we strongly encourage the usage of provided docker image as it would significantly reduce potential errors in the entire installation and setup process. (https://hub.docker.com/r/tjcadd2022/xmarkerfinder)
#### 3. Can I access and modify the codes used in xMarkerFinder?  
Yes. The codes used in xMarkerFinder are deposited in our GitHub repository and can be freely downloaded and modified according to users’ specific needs. However, the modification might cause unprecedented errors and we encourage users to try different parameters first, and then modify the codes. (https://github.com/tjcadd2020/xMarkerFinder)
#### 4. Can I use only certain steps of xMarkerFinder and skip other parts?  
Yes. The whole xMarkerFinder workflow contains four stages (23 steps) and every stage/step can be conducted independently and users could skip any one of them according to specific study designs.
#### 5. Can I use xMarkerFinder on environmental microbiome researches?  
Yes. Although xMarkerFinder is developed for human microbiome studies, it is also generalizable to other microbial habitats. 
#### 6. How long does it take to run xMarkerFinder?  
The time needed for the whole workflow depends on the dataset size, selected algorithm, and computational resources available. 
#### 7. What skills are required to run xMarkerFinder?  
Preliminary understanding of shell scripts would allow users to complete the whole workflow. Intermediate experience of R and python would facilitate users to better interpret and modify the codes.
#### 8. Is xMarkerFinder a pipeline for meta-analysis?  
Yes. xMarkerFinder aims to integrate different datasets and establish replicable markers. However, xMarkerFinder differs from systematic review as it integrates original datasets instead of the respective results.
### Part II Data processing
#### 1.	What kind of data should I use for xMarkerFinder?
Processed microbial count matrices and corresponding metadata are required. For cross-cohort analysis, we require merged datasets from at least three cohorts in the dicovery set to accomplish the full protocol with internal validations. xMarkerFinder is well adapted to microbial taxonomic and functional profiles derived from both amplicon and whole metagenomics sequencing data, as well as other omics layers, including but not limited to metatranscriptomics, metaproteomics, and metabolomics.
#### 2. If I don’t have the corresponding metadata, can I still use xMarkerFinder?
To perform meta-analysis, corresponding sample groups are required. Other metadata indices, such as body mass index, age and gender are recommended but not necessary.
#### 3.	Why should I normalize my data?
To mitigate challenges induced by different number of sequencing (e.g. library sizes), microbial count profiles are converted to relative abundances for subsequent analysis in xMarkerFinder.
#### 4.	Why should I perform data filtering?
To identify a replicable panel of microbial markers, we need to exclude rare microbial features, those with low occurrence rates across cohorts as they are not ideal candidates as global markers.
#### 5.	What does the training and test set do and why should I separate them?
To ensure models’ reliability, datasets are split to training and test set. Training set is used to train and have the model learn the hidden pattern. Test set is used to test the model after completing the training process and provides unbiased final model performance results.  
### Part III Using xMarkerFinder
#### 1.	How to solve installation errors?
Potential installation problems and solutions are provided along in our manuscript, and most problems would be avoided by simply using the docker image we provided instead of running all scripts locally (https://hub.docker.com/r/tjcadd2022/xmarkerfinder).
#### 2.	What machine learning classifier should I choose?
Step 6 provides the evaluation of multiple commonly used algorithms in machine learning, and users could choose the most suitable algorithm based on these results. However, due to its robustness and interpretability, Random Forest classifiers are considered suitable for most microbiome datasets. Therefore, step 6 is not compulsory and we recommend users to build Random Forest models first, move to other classifiers if they underperform.
#### 3.	How to choose suitable parameters when running xMarkerFinder?
For most scenarios, the default parameters would work. For further exploration, users are encouraged to try different parameters to get better results.
#### 4.	What is an AUC and how to interpret it?
AUC is the area under the ROC curve (the plot of the Sensitivity as y-axis versus 1-Specificity as x-axis). A perfect classifier gives an AUC of 1 while a simple classifier that makes completely random guesses gives an AUC of 0.5.

