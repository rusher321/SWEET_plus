# SWEEET
Sample-specific weighted correlation network (SWEET) method is desinged to model SINs by integrating the genome-wide sample weight with the differential correlation between the perturbed and aggregate networks.

## Input File Formats
- Gene expression profile (tab-delimited):
    * Column: Samples
    * Row: Genes
- Samples of interest: seperate with `\n`
- Genes of interest: seperate with `\n`

## Dependencies
The code is written in Python3. Additionally, the following package must also be installed:
- Numpy

## Basic Usage

Step 1: Calculate genome-wide sample weight:
```
python3 1.SWEET_sample_weight_calculating.py -f ./example/expression.txt -s ./example/weight.txt
```

`-h`: Get help for any of the commands
`-f`: Gene expression profile
`-k`: Balance parameter
`-s`: Output file name

Step 2: Calculate raw confidence scores of edges between given genes:
```
python3 2.SWEET_edge_weight_calculating.py -f ./example/expression.txt -w ./example/weight.txt -p ./example/patient.txt -g ./example/gene.txt -s ./example
```

`-h`: Get help for any of the commands
`-f`: Gene expression profile
`-w`: Sample weight file from previous step
`-p`: Patient file
`-g`: Gene file
`-s`: A path to output file

Step 3: Calculate the significance level of the confidence score for the edge between any two genes by a z-test:
```
python3 3.SWEET_calculating_mean_std.py -p ./example/patient.txt -l  ./example -s ./example/mean_std.txt -z False
```

`-h`: Get help for any of the commands
`-p`: Patient file
`-l`: Patient raw edge weight file path from previous step
`-s`: Output file name
`-z`: Calculate z-score or not

