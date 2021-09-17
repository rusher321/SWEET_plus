# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 13:16:09 2021

@author: sean
"""

import argparse
import numpy as np
import sys
import os

def check_file(expres) :
    checkset = set(["","NA","Na","na","nan","null"])
    for c in checkset :
        loc = np.where(expres==c)
        if loc[0].size :
            expres[loc] = "0"
            print(f"{c} in expres and transfer to 0")
    return expres

parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-p", type=str , default="./example/patient.txt" , help="patient file")
parser.add_argument("-l", type=str , default="./example" , help="file path")
parser.add_argument("-s", type=str , default="./example/mean_std.txt" , help="name of save file")

args = parser.parse_args()
file_p , file_l = args.p , (args.l).rstrip('/')
save = args.s

patlist = []
with open(file_p,mode='r') as rline :
    for nline in rline :
        tem = nline.strip('\n').split('\t')
        patlist.append(tem[0])

geneset = set()
pair = []
file = f"{file_l}/{patlist[0]}.txt"
if not os.path.exists(file) :
    print(f"{file} not found")
    sys.exit()
with open(f"{file_l}/{patlist[0]}.txt",mode='r') as rline :
    _ = rline.readline()
    for nline in rline :
        g1 , g2 , v , *_ = nline.strip('\n').split('\t')
        geneset.add(g1+'\t'+g2)
        pair.append(v)

for p in patlist[1:] :
    file = f"{file_l}/{p}.txt"
    if not os.path.exists(file) :
        print(f"{file} not found")
        sys.exit()
    with open(f"{file_l}/{p}.txt",mode='r') as rline :
        _ = rline.readline()
        for nline in rline :
            g1 , g2 , v , *_ = nline.strip('\n').split('\t')
            if (g1+'\t'+g2) not in  geneset :
                print(f"patient {p} gene pairs not map to other")
                pair.append(v)
pair = np.array(pair)
pair = check_file(pair)
pair = pair.astype(float)
vmean , vstd = np.mean(pair) , np.std(pair)

with open(save,mode='w') as wline :
    wline.write(f"mean\t{vmean}\nstd\t{vstd}\n")

print("Finish")