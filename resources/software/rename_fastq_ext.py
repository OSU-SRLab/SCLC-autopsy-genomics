#!/usr/bin/env python

#necessary since trim_galore outputs with _val_1.fq or _val_2.fq, and we want .fastq for consistency

import os
import re
import shutil
import sys

ext_re = re.compile(r"\.fastq$")
def rename_ext_in_path(orig_file_path, dir_with_file):
    #orig_file_path is properly named
    new_base = os.path.basename(orig_file_path)
    #either 1 or 2, from R1 or R2
    read_string = new_base.split("_")[3][1]
    new_path = os.path.join(dir_with_file, new_base)
    #derive trim_galore's output name
    base = ext_re.sub("_val_" + read_string + ".fq", new_base)
    path = os.path.join(dir_with_file, base)
    shutil.move(path, new_path)

if __name__ == "__main__":
    #from snakemake, the input file
    orig_fastq1_path = sys.argv[1]
    orig_fastq2_path = sys.argv[2]
    #the output path, i.e. where the files we actually want to rename are
    path_with_files = sys.argv[3]
    rename_ext_in_path(orig_fastq1_path, path_with_files)
    rename_ext_in_path(orig_fastq2_path, path_with_files)
