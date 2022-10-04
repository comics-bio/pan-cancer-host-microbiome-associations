# pan-cancer-host-microbiome-associations
Manuscript for pan-cancer-host-microbiome-associations


## Requirements

- This project requires R version __4.2.1__ (2022-06-23 ucrt) -- "Funny-Looking Kid"

- Please make sure that you have all required packages installed. You can
easily check whether that is the case by running

    ```bash
    conda create -n pancancer
    conda activate pancancer
    conda install -c conda-forge r-base=4.2.1
    Rscript requirements.R
    ```
    and confirm that it does not produce any errors.

## Workflow

To reproduce the main figures and result  of the project, you can follow these instructions.

### 0. Download data

#### RNA expression data matrix
The matrix `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` was generated following the Firehose pipeline: MapSplice + RSEM, then normalised by setting the upper-quartile to 1,000.

Pipeline details [here](https://gdc.cancer.gov/about-data/publications/pancanatlas) and here.

This was discussed in another thread here.

```bash
mkdir -p data/raw/PanCanAtlas/
cd data/raw/PanCanAtlas/
wget http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611
```
#### Microbiome data

```bash
mkdir -p data/raw/Kinght_2020
cd data/raw/Kinght_2020
wget 
```

### 1. Setup

```bash
cd ./src
Rscript preprocess_data.R
```

### SparseCCA

```bash

Options:
        -p PROJECT, --project=PROJECT
                TCGA cancer project [default: TCGA-BRCA]

        -o OUTPUT, --output=OUTPUT
                The analysis result location [default: ../result/]

        -n NUM, --num=NUM
                Cluster cores numbers [default: 12]

        -c COMPONET, --componet=COMPONET
                CCA componet[default: 10]

        -f FDR, --fdr=FDR
                CCA componet[default: 0]

        -h, --help
                Show this help message and exit
```


### Lasso 


This script implements a lasso regression for each host gene’s expression as response and abundances of microbial taxa and values of other covariates as independent variables. It uses leave-one-out cross-validation to estimate the tuning parameter, which is used to fit the final model on a given disease dataset. It uses desparsified lasso approach (R package HDI) to obtain 95% confidence intervals and p-values for the coefficient of each microbe associated with a given host gene. This script implements a parallel framework for executing the gene-wise lasso analysis, where execution of lasso models on host genes can be parallelized across multiple nodes and cores on a compute cluster.




## Contact

If you have questions about the code or the analyses presented in this project,
please feel free to contact me, 12233060@mail.sustech.edu.cn
:smiley:
