# SWEEET
Here, we propose a sample-specific weighted correlation network (SWEET) method to model SINs by integrating the genome-wide sample weight with the differential correlation between the perturbed and aggregate networks.

## Input File Formats
- gene expression profile:
    * column: patients(samples)
    * row: genes
- interest patients(samples): seperate with `\n`
- interest genes: seperate with `\n`

## Basic Usage

To run step 1 for calculte sample weight:
```
python3 1.SWEET_sample_weight_calculating.py -f ./example/expression.txt -s ./example/weight.txt
```

`-h`: Get help for any of the commands
`-f`: Gene expression profile
`-k`: Balance parameter
`-s`: Output file name

To calculate raw edge score between genes:
```
python3 2.SWEET_edge_weight_calculating.py -f ./example/expression.txt -w ./example/weight.txt -p ./example/patient.txt -g ./example/gene.txt -s ./example
```

`-h`: Get help for any of the commands
`-f`: Gene expression profile
`-w`: Sample weight file from previous step
`-p`: Patient file
`-g`: Gene file
`-s`: A path to output file

To calculate the statical significance:
```
python3 3.SWEET_calculating_mean_std.py -p ./example/patient.txt -l  ./example -s ./example/mean_std.txt -z False
```

`-h`: Get help for any of the commands
`-p`: Patient file
`-l`: Patient raw edge weight file path from previous step
`-s`: Output file name
`-z`: Calculate z-score or not