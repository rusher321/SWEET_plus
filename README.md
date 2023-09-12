# SWEET_plus
Sample-specific weighted correlation network (SWEET) method is desinged to model SINs by integrating the genome-wide sample weight with the differential correlation between the perturbed and aggregate networks.
SWEET_plus is an expanded version of SWEET for modeling SINs of two omics-data, eg. bacteria-virome or others. 

## Input File Formats
- Gene expression matrix (tab-delimited):
    * Column: Samples
    * Row: Genes
- Samples of interest: seperate with `\n`
- Genes of interest: seperate with `\n`

## Dependencies
The code is written in Python3. Additionally, the following package must also be installed:
- Numpy
- Pandas
- scipy
- argparse


## Basic Usage
The example datasets are stored inside example folder, as well the example outputs.  
Step 1: calculate genome-wide sample weight:

*** raw method ***
```
python3 1.SWEET_sample_weight_calculating.py -f ./example/expression.txt -s ./example/weight.txt
```

*** use two matrix as the input and spearman correlation ***
```
python3 1.SWEET_sample_weight_calculating.py -f ./example/expression.txt -f2 ./example/expression2.txt -s ./test/weight.spearman -c spearman
```
  
  `-h`: --help  show this help message and exit
  `-f`: F        A path to 'gene expression matrix' or 'gene/species correlation matrix' file
  `-f2`: F2      A path to 'gene expression matrix2/species correlation matrix' file
  `-k`: K        balance parameter
  `-n`: N        network or expession matrix
  `-c`: C        correlation method [spearman/pearson]
  `-s`: S        A path to the output 'sample weight' file

Step 2: calculate confidence scores of edges between given genes for each sample of interest:

*** raw method ***
```
python3 2.SWEET_edge_score_calculating.py -f ./example/expression.txt -w ./example/weight.txt -p ./example/patient.txt -g ./example/gene.txt -s ./example
```
*** use two matrix as the input and spearman correlation ***
```
python3 2.SWEET_edge_score_calculating.py -f ./example/expression.txt -f2 ./example/expression2.txt -w ./test/weight.spearman.mean.txt -p ./example/patient.txt -g ./example/gene.txt -c spearman -s ./test/
```

  `-h`: --help  show this help message and exit
  `-f`: F        A path to 'gene expression matrix' file
  `-f2` F2      A path to 'gene expression matrix2/species correlation matrix' file
  `-w` W        A path to 'sample weight' file (i.e., the output file from step1)
  `-p` P        A path to 'samples of interest' file
  `-g` G        A path to 'genes of interest' file
  `-c` C        correlation method [spearman/pearson]
  `-s` S        A path to the output 'confidence scores of edges' files for each sample of interest

Step 3: calculate the significance level of the confidence score for the edge between any two genes by a z-test:
```
python3 3.SWEET_calculating_mean_std_zscore.py -p ./example/patient.txt -l  ./example -s ./example/mean_std.txt -z False
```

`-h`: Get help with the commands  
`-p`: A path to "samples of interest" file  
`-l`: A path to the "confidence scores of edges" file for each sample of interest (i.e., the output files from step 2)  
`-s`: A path to the output file(s)  
`-z`: Indicates whether the calculation of z score (Ture) or not (False)  

Note that the mean and standard deviation are calculated by the confidence scores of all edges for the samples of interest; therefore, different lists of "samples of interest" will generate distinct means and standard deviations.

