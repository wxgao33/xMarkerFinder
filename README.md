# Identification and validation of microbial marker from cross-cohort datasets using xMarkerFinder
xMarkerFinder is a four-stage workflow for microbiome research including cross-cohort differential signature analysis, model construction, model validation, and marker identification.  



## User Tutorial
#### 1
  `$ Rscript 1_Normalization.R -W /workplace/ -p abundance.txt -o TEST`

## FAQs
### Part I General questions
#### 1. When should I use xMarkerFinder?  
xMarkerFinder is suitable for microbial marker identification from cross-cohort datasets. Our previous studies demonstrated its applicability in identifying global microbial diagnostic markers for adenoma and colorectal cancer. Moreover, xMarkerFinder could also be applied to marker determination in disease prognosis, treatment stratification, metastasis surveillance, adverse reactions anticipation etc. Any research dedicated to marker identification from multi-population microbial datasets are welcome.
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
Processed microbial absolute abundance matrices and corresponding metadata are required. For cross-cohort analysis, we require merged datasets from at least three cohorts to accomplish the full protocol with internal validations. xMarkerFinder is well adapted to microbial taxonomic and functional profiles derived from both amplicon and whole metagenomics sequencing data, as well as other omics layers, including but not limited to metatranscriptomics, metaproteomics, and metabolomics.
#### 2. If I don’t have the corresponding metadata, can I still use xMarkerFinder?
To perform meta-analysis, corresponding sample groups are required. Other metadata indices, such as body mass index, age and gender are recommended but not necessary.
#### 3.	Why should I perform data filtering?
To identify a replicable panel of microbial markers, we need to exclude rare microbial features, those with low occurrence rates and relative abundances across cohorts as they are not ideal candidates as global markers.
#### 4.	Why should I normalize my data?
To mitigate challenges induced by different number of sequencing (e.g. library sizes), microbial profiles are converted to relative abundances for subsequent analysis in xMarkerFinder.
#### 5.	Why should I perform data filtering?
To identify a replicable panel of microbial markers, we need to exclude rare microbial features, those with low occurrence rates and relative abundances across cohorts as they are not ideal candidates as global markers.
#### 6.	What does the training and test set do and why should I separate them?
To ensure models’ reliability, datasets are split to training and test set. Training set is used to train and have the model learn the hidden pattern. Test set is used to test the model after completing the training process and provides unbiased final model performance metric results.  
### Part III Using xMarkerFinder
#### 1.	How to solve installation errors?
Potential installation problems and solutions are provided along in our manuscript, and most problems would be avoided by simply using the docker image we provided instead of running all scripts locally.
#### 2.	What machine learning classifier should I choose?
Step 6 provides the evaluation of multiple commonly used algorithms in machine learning, and users could choose the most suitable algorithm based on these results. However, due to its robustness and interpretability, Random Forest classifiers are deemed suitable for most microbiome datasets. Therefore, step 6 is not compulsory and we recommend users to build Random Forest models first, move to other classifiers if they do not perform well.
#### 3.	How to choose suitable parameters when running xMarkerFinder?
For most scenarios, the default parameters would work. For further exploration, users are encouraged to try different parameters to get better results.
#### 4.	What is an AUC and how to interpret it?
AUC is the area under the ROC curve (the plot of the Sensitivity as y-axis versus 1-Specificity as x-axis). A perfect classifier gives an AUC of 1 while a simple classifier that makes completely random guesses gives an AUC of 0.5.

