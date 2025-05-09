#!/usr/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import re
import subprocess

# /bi/apps/TrAELseq/latest/TrAEL-seq/rename_trael_samples.py --index 135 --sample_names "D5, D20, D40"
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description = '''For adding in sample names to TrAEL files. this will be applied ot all files in the current working directory''')
parser.add_argument('--sample_names', type=str, default="", help='Comma separated set of names - must be in correct order to match the supplied indices')
parser.add_argument('--index', type=str, default="", help='Indexes to rename')

args=parser.parse_args()

index_dict = {
    "1": "_AGTC_index1_",
    "2": "_GACT_index2_",
    "3": "_CTTG_index3_",
    "4": "_TCGA_index4_",
    "5": "_AAGG_index5_",
    "6": "_TTCC_index6_",
    "7": "_GTGC_index7_",
    "8": "_GCCA_index8_",
    "9": "_GATG_index9_"
}

# create sample name dictionary
sample_dict = {}

def main():

    indexes = list(args.index.strip())

    sample_names = args.sample_names.strip().split(",")

    if (len(indexes) == len(sample_names)):

        count = 0

        for index in indexes:
            #print(f"sample name = {sample_names[count].strip()}")

            cleaned_name = clean_sample_name(sample_names[count].strip())
            sample_dict[index] = cleaned_name
            #print(f"cleaned name = {cleaned_name}")

            old_name = index_dict[index]
            new_name = f"{old_name}{cleaned_name}_"
            cmd = f"rename {old_name} {new_name} *{old_name}*"
            print(cmd)

            try:
                subprocess.run(cmd, shell=True, executable="/bin/bash")

            except Exception as err:
                print(f"\n !! Couldn't run rename command !!")
                print(err)

            count+=1

    else: 
        print("!! number of indexes and sample names do not match - please check and try again!!")
        print(f"number of sample names: {len(sample_names)}, number of indexes: {len(indexes)}")
        exit()


#---------------------------------------------------
# remove any unwanted characters from sample names
#---------------------------------------------------
def clean_sample_name(sample_name):    
    return re.sub(r'[^a-zA-Z0-9.\-_]+_?', '_', sample_name)


if __name__ == "__main__":
    main()