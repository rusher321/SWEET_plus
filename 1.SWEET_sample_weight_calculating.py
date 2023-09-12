import os
import argparse
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-f", type=str , default="./example/expression.txt" , help="A path to 'gene expression matrix' or 'gene/species correlation matrix' file")
parser.add_argument("-f2", type=str , default="./example/expression2.txt" , help="A path to 'gene expression matrix2/species correlation matrix' file", required = False)
parser.add_argument("-k", type=float , default=0.1 , help="balance parameter")
parser.add_argument("-n", type = bool, default = False, help="network or expession matrix")
parser.add_argument("-c", type = str, default = "pearson", help = "correlation method [spearman/pearson]")
parser.add_argument("-s", type=str , default="./example/weight.txt" , help="A path to the output 'sample weight' file")

args = parser.parse_args()
#global pat
#global patlen

def read_file(file):
    gene, value = [], []
    with open(file,mode='r') as rline :
        global pat
        pat=rline.readline().strip('\n').split('\t')[1:]
        global patlen
        patlen=len(pat)
        for nline in rline :
            g , *v = nline.strip('\n').split('\t')
            value += v
            gene.append(g)
    genelen = len(gene)
    value = np.array(value,dtype=float).reshape(genelen,patlen).T
    return value

def weight_comput(file, title):
    value = read_file(file)
    if args.n:
        value = (np.sum(value,axis=1)-1)/(patlen-1)
        rmax , rmin = np.max(value) , np.min(value)
        dif = rmax - rmin + 0.01
        value = (value - rmin + 0.01)/dif
        value = value * k * patlen
    else:
        if(args.c == "pearson"):
            value = np.corrcoef(value)
        elif(args.c == "spearman"):
            value = stats.spearmanr(value.T)
            value = value.correlation
        else:
            print("correlation method need to be use pearson or spearman!")
            sys.exit()
        value = (np.sum(value,axis=1)-1)/(patlen-1)
        rmax , rmin = np.max(value) , np.min(value)
        dif = rmax - rmin + 0.01
        value = (value - rmin + 0.01)/dif
        value = value * k * patlen
    with open(save+title+".txt", mode='w') as w_line :
            w_line.write('patient\tsample_weight\n')
            for p,v in zip(pat,value) :
                tem = p + '\t' + str(v) + '\n'
                w_line.write(tem)
    return value

def main(args):
    f1 = weight_comput(args.f, title = ".matrix1")
    if os.path.exists(args.f2):
        f2 = weight_comput(args.f2, title = ".matrix2")
        f3 = f1+f2
        with open(save+".mean.txt", mode='w') as w_line :
            w_line.write('patient\tsample_weight\n')
            for p,v in zip(pat, f3) :
                tem = p + '\t' + str(v) + '\n'
                w_line.write(tem)

k, save = args.k, args.s

if __name__ == "__main__":
    main(args)
 
print("Finish")
