# Microbial marker identification from cross-cohort datasets (MMIC)
MMIC is a curated four-stage workflow for microbiome researches including cross-cohort differential signature analysis, model construction, model validation and marker identification.

## FAQ
### Part I General questions
#### 1. When should I use MMIC?  
MMIC is suitable for microbial marker identification from cross-cohort datasets. Our previous studies demonstrated its applicability in identifying global microbial diagnostic markers for adenoma and colorectal cancer. Moreover, MMIC could also be applied to marker determination in disease prognosis, treatment stratification, metastasis surveillance, adverse reactions anticipation etc. Any research dedicated to marker identification from multi-population microbial datasets are welcome.
#### 2. How should I setup the required computational environment for MMIC?  
We provide detailed instructions on software installation for users to run the whole MMIC workflow locally. However, we strongly encourage the usage of provided docker image as it would significantly reduce potential errors in the entire installation and setup process. (https://hub.docker.com/r/tjcadd2022/mmic)
#### 3. Can I access and modify the codes used in MMIC?  
Yes. The codes used in MMIC are deposited in our GitHub repository and can be freely downloaded and modified according to users’ specific needs. However, the modification might cause unprecedented errors and we encourage users to try different parameters first, and then modify the codes. (https://github.com/tjcadd2020/MMIC)
#### 4. Can I use only certain steps of MMIC and skip other parts?  
Yes. The whole MMIC workflow contains four stages (23 steps) and every stage/step can be conducted independently and users could skip any one of them according to specific study designs.
#### 5. Can I use MMIC on environmental microbiome researches?  
Yes. Although MMIC is developed for human microbiome studies, it is also generalizable to other microbial habitats. 
#### 6. How long does it take to run MMIC?  
The time needed for the whole workflow depends on the dataset size, selected algorithm, and computational resources available. 
#### 7. What skills are required to run MMIC?  
Preliminary understanding of shell scripts would allow users to complete the whole workflow. Intermediate experience of R and python would facilitate users to better interpret and modify the codes.
#### 8. Is MMIC a pipeline for meta-analysis?  
Yes. MMIC aims to integrate different datasets and establish replicable markers. However, MMIC differs from systematic review as it integrates original datasets instead of the respective results.


