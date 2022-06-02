# Identification and validation of microbial marker from cross-cohort datasets using xMarkerFinder
xMarkerFinder is a four-stage workflow for microbiome research including cross-cohort differential signature analysis, model construction, model validation, and marker identification. Detailed [scripts](./scripts), [example files](./data), and a ready-to-use [docker image](https://hub.docker.com/repository/docker/tjcadd2022/xmarkerfinder) are provided.  

![ ](https://img.shields.io/badge/python-3.7-blue) ![GitHub top language](https://img.shields.io/github/languages/top/tjcadd2020/xMarkerFinder)  ![GitHub](https://img.shields.io/github/license/tjcadd2020/xMarkerFinder) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/tjcadd2020/xMarkerFinder) [![GitHub version](https://badge.fury.io/gh/tjcadd2020%2FxMarkerFinder.svg)](https://badge.fury.io/gh/tjcadd2020%2FxMarkerFinder) ![GitHub issues](https://img.shields.io/github/issues/tjcadd2020/xMarkerFinder) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tjcadd2020/xMarkerFinder/HEAD) [![](https://img.shields.io/badge/website-CADD-lightgrey)](https://cadd.tongji.edu.cn/)

## Table of Contents
* [User Tutorial](#user-tutorial)
  * [Stage 1 Cross-cohort differential signature analysis](#stage-1-cross-cohort-differential-signature-analysis)
  * [Stage 2 Model construction](#stage-2-model-construction)
  * [Stage 3 Model validation](#stage-3-model-validation)
  * [Stage 4 Marker interpretation](#stage-4-marker-interpretation)
* [FAQs](#faqs)
  * [Part I General questions](#part-i-general-questions)
  * [Part II Data processing](#part-ii-data-processing)
  * [Part III Using xMarkerFinder](#part-iii-using-xmarkerfinder)


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
### Stage 3 Model validation
#### 13.	Internal validations (intra-cohort, cohort-to-cohort, and LOCO validation). 
As stated above, this step provides extensive internal validations to ensure the robustness and reproducibility of identified markers in different cohorts via intra-cohort validation, cohort-to-cohort transfer, and LOCO validation. Output files contain multiple performance metrics used to assess the markers internally, including AUC, specificity, sensitivity, accuracy, precision and F1 score.  
```
$ python 13_Validation.py -W /workplace/ -m metadata.txt -p candidate_marker.txt -g Group -e exposure -b Cohort -c classifier -s 0 -o TEST
```
```
-p input optimal candidate marker file (output file of Step 9 or Step 10)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers.  
- Output files:  
validation.txt: the overall performance of markers in internal validations.  
#### 14.	Internal Validation grid plot. 
This step provides the visualization of Step 13. Users could select any output file of Step 13 and draw a grid plot.  
```
$ python 14_Validation_grid_plot.py -W /workplace/ -p validation_metric.txt -o TEST
```
```
-p input tableau file for visualization (output file of Step 13)
```
- Input files:
validation_metric.txt: the overall performance of markers in internal validations.
- Output files:
validation.pdf: the visualization of input file.
#### 15.	External test. 
As the best-performing candidate markers and classification model are established, the test dataset is used to externally validate their generalizability. The input external metadata and microbial relative profiles need to be in the same format as initial input files for the training dataset. This step returns the overall performance of the model and its AUC plot.  
```
$ python 15_Test.py -W /workplace/ -m metadata.txt -p candidate_marker.txt -a external_metadata.txt -x external_profile.txt -g Group -e exposure -c classifier -r best_param.txt -s 0 -o TEST
```
```
-a input external metadata file for the test dataset
-x input external microbial relative abundance file as the test dataset
-r input optimal hyperparameter file (output file of Step 11)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers.  
external_metadata.txt: the clinical metadata of the external test dataset.  
external_profile.txt: the relative abundance matrix of the external test dataset. 
best_param.txt: the best hyperparameter combination of classification model.  
- Output files:  
test_result.txt: the overall performance of model in external test dataset.  
test_auc.pdf: the visualization of the AUC value in test_result.txt.  
#### 16.	Specificity assessment.   
To further assess markers’ specificity for experimental group of interest, they are used to construct classification models to discriminate between other related diseases and corresponding controls. Cross-validation AUC values of other classification models and visualization are returned.   
```
$ python 16_Specificity.py -W /workplace/ -p candidate_marker.txt -a other_metadata.txt -x other_profile.txt -g Group -e exposure -b Cohort -c classifier -r best_param.txt -s 0 -o TEST
```
```
-a input metadata file of samples from other diseases
-x input microbial relative abundance file of samples from other diseases
-e the control group name (in example file: CTR)
-b the column name of cohort(in example file: Cohort)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers.  
other_metadata.txt: the clinical metadata of samples for other diseases.  
other_profile.txt: the relative abundance matrix of other diseases.  
- Output files:  
specificity_result.txt: AUC values of models constructed with candidate markers in other related diseases.  
specificity_auc.pdf: the visualization of the specificity_result.txt.  
#### 17.	Alternative specificity assessment.   
Random samples of case and control class of other diseases are added into the classification model, respectively, both labelled as “control”, the variations of corresponding AUCs of which are calculated used for visualization.   
```
$ python 17_Specificity_add.py -W /workplace/ -m metadata.txt -p candidate_marker.txt -q external_metadata.txt -l external_profile.txt -a other_metadata.txt -x other_profile.txt -g Group -e exposure -b Cohort -c classifier -r hyperparamter.txt -n 5 -s 0 -o TEST
```
```
-q input external metadata file for the test dataset
-l input external microbial relative abundance file as the test dataset
-a input metadata file of samples from other diseases
-x input microbial relative abundance file of samples from other diseases
-e the control group name (in example file: CTR)
-b the column name of cohort(dataset)
-n the number of samples to add into the model each time 
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers.  
external_metadata.txt: the clinical metadata of the external test dataset.  
external_profile.txt: the relative abundance matrix of the external test dataset.  
other_metadata.txt: the clinical metadata of samples for other diseases.  
other_profile.txt: the relative abundance matrix of other diseases.  
- Output files:  
specificity_add_result.txt: AUC values of models constructed with candidate markers in other related diseases.  
specificity_add_auc.pdf: the visualization of the specificity_result.txt.  
### Stage 4 Marker interpretation.
#### 18.	Marker importance.
Permutation feature importance is employed here to evaluate markers’ contributions in the best-performing classification model.  
```
$ python 18_Marker_importance.py -W /workplace/ -m metadata.txt -p candidate_marker.txt -g Group -e exposure -c classifier -r best_param.txt -s 0 -o TEST
```
```
-p input candidate markers (output file of Step 9 or Step 10)
-r input optimal hyperparameter file (output file of Step 11)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
candidate_marker.txt: the optimal panel of candidate markers.  
best_param.txt: the best hyperparameter combination of classification model.  
- Output files:  
marker_importance.txt: permutation feature importance of candidate markers via ten permutations.  
#### 19.	Marker importance plot.   
This step provides the visualization of Step 18.  
```
$ python 19_Marker_importance_plot.py -W /workplace/ -p marker_importance.txt -o TEST
```
```
-p input file for plot (output file of Step 18)
```
- Input files:  
metadata.txt: the clinical metadata of the training dataset.  
marker_importance.txt: permutation feature importance of candidate markers.  
- Output files:  
marker_importance.pdf: the visualization of feature importance file.  
#### 20.	Microbial taxon correlation.   
Inter-microbiota correlation is calculated using FastSpar with 50 iterations and the output files contain the correlation and P value between every microbiota pair.   
As the input file for Step 20 needs to be microbial count profile in .tsv format where each row describes a microbial taxon and each column represents a sample (could be converted profiles of all features, differential signatures, or candidate markers according to users’ need, and null values needed to be set as 0) and header needs to start with “#OTU ID”, an additional file conversion script is provided.  
```
$ python 20_Convert.py -W /workplace/ -p abundance.txt -s selected_feature.txt -o TEST
```
```
-p input feature raw count file before normalization.
-s selected features for calculating microbial correlation (could be differential signatures or candidate markers, output file of Step 4,9, or 10).
```
- Input files:  
abundance.txt: microbial raw count profile before normalization.  
selected_feature.txt: selected features for calculating microbial correlation (output file of Step 4, 9, or 10)  
- Output files:  
convert.tsv: the converted file appropriate for calculating microbial taxon correlation.  
```
$ ./20_Microbial_taxon_correlation.sh -W /workplace/ –i feature_abundance.tsv -o TEST -t 4
```
```
-i input feature abundance file  
-t threads of available computational source  
```
- Input files:  
microbial_taxon.tsv: microbial count profile in .tsv format where each row describes a microbial taxon and each column represents a sample and the header needs to start with “#OTU ID”. Example input file is provided and users are recommended to user the conversion script above to convert files into appropriate formats.  
-t the threads of available computational source when running  
- Output files:  
median_correlation.tsv: the correlation coefficients between every input taxon pair.  
pvalues.tsv: the statistical significance of median_correlation.tsv.  
#### 21. Microbial taxon correlation plot. 
The visualization of Step 20 is performed using Gephi.  
(i) Preprocess of the results of Step 20 to ensure that Step 21 only draws significant correlations (pvalues<0.05) with absolute correlation coefficients above 0.5 (default).
```
$ python 21_Microbial_taxon_correlation_plot.py -W /workplace/ –c median_correlation.tsv -p pvalues.tsv -t 0.5 -o TEST 
```
```
-c input correlation profile (output file of Step 20)
-p input pvalue profile (output file of Step 20)
-t input correlation threshold (default: 0.5)
```
- Input files:
median_correlation.tsv: the correlation coefficients profile (output file of Step 20).
pvalues.tsv: the statistical significance of median_correlation.tsv (output file of Step 20).
- Output files:
correlation.csv: adjusted correlation profile for Gephi input, only significant correlations reserved.  
(ii)	Open Gephi and click “File” – “Import spreadsheet”.
<img width="415" alt="image" src="https://user-images.githubusercontent.com/54845977/171319514-fe15e6d3-6d44-43c5-9c50-3e6c5d52a5e5.png">  
(iii) Choose the adjusted correlation profile.  
<img width="415" alt="image" src="https://user-images.githubusercontent.com/54845977/171319788-af5776f9-7b91-4dc8-a55f-d6ccc404a2c4.png">  
(iv)  Choose a preferable layout type to form the basic network and press the “stop” button when the network becomes stable (Fruchterman Reingold style is recommended).  
<img width="415" alt="image" src="https://user-images.githubusercontent.com/54845977/171319813-0fed579e-6c7d-4581-bf7e-174aa8d391e1.png">   
(v)  For further optimization of the network, appearances of nodes and edges should be adjusted according to users’ need, as well as the labels of nodes. <img width="415" alt="image" src="https://user-images.githubusercontent.com/54845977/171319835-a572168e-fad4-47d0-a03b-16b528c99d54.png">    

#### 22.	Multi-omics correlation. 
If users have multi-omics or multidimensional microbial profiles of the same dataset, the correlation between different omics or dimensions are calculated via HAllA.
```
$ ./22_Multi_omics_correlation.sh -W /workplace/ -i microbial_abundance_1.txt -d microbial_abundance_2.txt -o TEST
```
```
-i input microbial abundance file 1
-d input microbial abundance file 2
```
- Input files:  
microbial_abundance_1.txt: microbial abundance profile 1.  
microbial_abundance_2.txt: microbial abundance profile 2. These two input files should have the same samples (columns) but different features (rows).  
- Output files:  
results/all_associations.txt: associations between different omics or dimensions.  
#### 23.	Multi-omics correlation plot. 
The visualization step of Step 22.  
```
$ ./23_Multi_omics_correlation_plot.sh -W /workplace/ -i TEST
```
```
-i input file for plot (output file of Step 22)
```
- Input files:  
results: the result folder of Step 22  
- Output files:  
results/hallagram.png: the visualization of all_associations.txt with only significant associations highlighted.  











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
The time needed for the whole workflow depends on the dataset size, selected algorithm, and computational resources available. The following time estimates are based on execution of our protocol on provided example datasets with RF (Timing 1), DT (Timing 2), KNN (Timing 3), and SVC (Timing 4) using the xMarkerFinder docker image on a MacBook Pro (2.4-GHz quad-core eighth-generation Intel Core i5 processor, 16-GB 2133-MHz LPDDR3 memory).  
|     Step    |     Timing 1       |     Timing 2      |     Timing 3      |     Timing 4      |
|-------------|--------------------|-------------------|-------------------|-------------------|
|     1       |     0m20.600s      |     0m20.600s     |     0m20.600s     |     0m20.600s     |
|     2       |     0m11.372s      |     0m11.372s     |     0m11.372s     |     0m11.372s     |
|     3       |     1m21.356s      |     1m21.356s     |     1m21.356s     |     1m21.356s     |
|     4       |     0m22.535s      |     0m22.535s     |     0m22.535s     |     0m22.535s     |
|     5       |     0m2.323s       |     0m2.323s      |     0m2.323s      |     0m2.323s      |
|     6       |     0m12.464s      |     0m12.464s     |     0m12.464s     |     0m12.464s     |
|     7       |     1m43.332s      |     0m3.252s      |     0m3.105s      |     0m50.913s     |
|     8       |     0m0.863s       |     0m1.015s      |     0m1.178s      |     0m1.102s      |
|     9       |     18m37.552s     |     0m53.413s     |     0m21.103s     |     10m32.261s    |
|     10      |     24m59.785s     |     24m59.785s    |     24m59.785s    |     24m59.785s    |
|     11      |     8m42.417s      |     0m31.801s     |     0m39.348s     |     0m32.112s     |
|     12      |     0m15.083s      |     0m3.149s      |     0m3.564s      |     0m3.724s      |
|     13      |     123m20.805s    |     4m27.044s     |     6m6.515s      |     10m12.050s    |
|     14      |     0m5.467s       |     0m4.541s      |     0m4.655s      |     0m3.428s      |
|     15      |     0m7.120s       |     0m4.002       |     0m3.761s      |     0m3.150s      |
|     16      |     2m25.064s      |     0m5.315s      |     0m5.287s      |     0m4.426s      |
|     17      |     2m21.262s      |     0m6.646s      |     0m7.342s      |     0m22.211s     |
|     18      |     0m43.265s      |     0m2.041s      |     0m18.809s     |     0m14.746s     |
|     19      |     0m3.794s       |     0m2.938s      |     0m3.062s      |     0m2.264s      |
|     20      |     6m32.696s      |     6m32.696s     |     6m32.696s     |     6m32.696s     |
|     21      |     0m0.884s       |     0m0.884s      |     0m0.884s      |     0m0.884s      |
|     22      |     7m57.119s      |     7m57.119s     |     7m57.119s     |     7m57.119s     |
|     23      |     0m18.923s      |     0m18.923s     |     0m18.923s     |     0m18.923s     |
#### 7. What skills are required to run xMarkerFinder?  
Preliminary understanding of shell scripts would allow users to complete the whole workflow. Intermediate experience of R and python would facilitate users to better interpret and modify the codes.
#### 8. Is xMarkerFinder a pipeline for meta-analysis?  
Yes. xMarkerFinder aims to integrate different datasets and establish replicable markers. However, xMarkerFinder differs from systematic review as it integrates original datasets instead of the respective results.
### Part II Data processing
#### 1.	What kind of data should I use for xMarkerFinder?
Processed microbial count matrices and corresponding metadata are required. For cross-cohort analysis, we require merged datasets from at least three cohorts in the dicovery set to accomplish the full protocol with internal validations. xMarkerFinder is well adapted to microbial taxonomic and functional profiles derived from both amplicon and whole metagenomics sequencing data, as well as other omics layers, including but not limited to metatranscriptomics, metaproteomics, and metabolomics.
#### 2. If I don’t have the corresponding metadata, can I still use xMarkerFinder?
To perform meta-analysis, corresponding sample groups are required. Other metadata indices, such as body mass index, age and gender are recommended but not necessary. However, it is worth noticing that the absence of metadata information might compromise the correction for confounding effects and the identification of microbial markers.  
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

