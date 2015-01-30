#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
from collections import deque
from subprocess import call
from optparse import OptionParser
from tempfile import mkstemp
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import time
import resource
#import file


FNULL = open('/dev/null', 'w')
base_path = os.path.dirname(sys.argv[0])[:-len('src/')]


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


def compress(options):
    """
    Compress the reads using all methods.
    """

    ensure_dir(options.output_dir + '/goodbad/')
    ensure_dir(options.output_dir + '/original/')

    std_err_file = open('compress.log', 'w')

    GB_COMPRESSION_CMD = "./src/good_bad_coding.py -r [READ] -c 2 -b 0 "
    
    # Store which compression directories we created.
    options.compressed_dirs = []
    options.compressed_dirs.append('original')
    options.compressed_dirs.append('goodbad')

    for reads_filename in options.reads_filenames:
        
        # Copy the original sequences over.
        out_cmd("", std_err_file.name, ["cp", reads_filename, options.output_dir + '/original/' + os.path.basename(reads_filename)])
        shutil.copyfile(reads_filename, options.output_dir + '/original/' + os.path.basename(reads_filename))

        # Good/bad binary compression.
        call_arr = GB_COMPRESSION_CMD.replace('[READ]', reads_filename).split()
        output_fp = open(options.output_dir + '/goodbad/' + os.path.basename(reads_filename), 'w')

        out_cmd(options.output_dir + '/goodbad/' + os.path.basename(reads_filename), std_err_file.name, call_arr)
        call(call_arr, stdout=output_fp, stderr=std_err_file)

        # Polynomial regression.

        # Profile regression.

    std_err_file.close()

    pass


def decompress(options):
    """
    After compressing the reads, decompress them. so they can be used in downstream analyses.
    """
    pass


def calc_mean_squared_error(options):
    """
    Calculate mean squared error between the original and decompressed reads.
    """
    pass


def quality_preprocessing(options):
    """
    Examine the effects of lossy compression on quality preprocessing tools.
    """
    pass


def assemble(options):
    """
    Test assemblies using ALLPATHS-LG.
    """

    std_err_file = open('compress.log', 'a')

    # The first thing step is to create the in_groups.csv and in_libs.csv.

    IN_GROUPS_CSV = """group_name, library_name,   file_name
frag,   Illumina_01,   [FULL_PATH]/[COMPRESSION]/frag_*.fastq 
shortjump,  Illumina_02,    [FULL_PATH]/[COMPRESSION]/shortjump_*.fastq"""

    IN_LIBS_CSV = """library_name,   project_name,   organism_name,  type,   paired, frag_size,  frag_stddev,    insert_size,    insert_stddev,  read_orientation,   genomic_start,  genomic_end
Illumina_01,    assembly,    unknown,    fragment,   1,  180,    10, ,   ,   inward, ,   
Illumina_02,    assembly,    unknown,    jumping,    1,  ,   ,   3000,   500,    outward,    ,"""

    #print(IN_GROUPS_CSV.replace('[FULL_PATH]', os.path.abspath(options.output_dir)).replace('[COMPRESSION]', 'goodbad'))

    for compression_method in options.compressed_dirs:
        ensure_dir(options.output_dir + '/assemble/' + compression_method + '/')

        open(options.output_dir + '/assemble/' + compression_method + '/in_groups.csv', 'w').write(IN_GROUPS_CSV.replace('[FULL_PATH]', \
                os.path.abspath(options.output_dir)).replace('[COMPRESSION]', compression_method))
        open(options.output_dir + '/assemble/' + compression_method + '/in_libs.csv', 'w').write(IN_GROUPS_CSV.replace('[FULL_PATH]', \
                os.path.abspath(options.output_dir)).replace('[COMPRESSION]', compression_method))

        ALLPATHS_CMD = "RunAllPathsLG PRE=" + options.output_dir + '/assemble/' + compression_method + " DATA_SUBDIR=. RUN=allpaths SUBDIR=run THREADS=32 OVERWRITE=True"
        # RunAllPathsLG PRE=. REFERENCE_NAME=. DATA_SUBDIR=. RUN=allpaths SUBDIR=run THREADS=32 OVERWRITE=True
        # RunAllPathsLG PRE=/assemblies DATA=datadir RUN=allpaths SUBDIR=run THREADS=32 OVERWRITE=True

        call_arr = ALLPATHS_CMD.split()
        #output_fp = open(options.output_dir + '/goodbad/' + os.path.basename(reads_filename), 'w')

        out_cmd("", std_err_file.name, call_arr)


def align_reads(options):
    """

    """
    pass


"""
I/O Helpers
"""
def setup_shell_file():
    if shell_file_fp:
        shell_file_fp.write("#!/bin/bash\n")

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
    assert os.path.exists(d)

def out_cmd(std_out = "", std_err = "", *objs):
    #line(75)
    if shell_file_fp:
        if std_out != "":
            std_out_sht = " 1>%s " % (std_out)
        else:
            std_out_sht = ""
        if std_err != "":
            std_err_sht = " 2>%s " % (std_err)
        else:
            std_err_sht = ""
        shell_file_fp.write(' '.join(*objs) + std_out_sht + std_err_sht + "\n")
        shell_file_fp.flush()
    print(bcolors.OKBLUE + "COMMAND:\t" + bcolors.ENDC, ' '.join(*objs) + std_out_sht, file=sys.stderr)
"""
I/O Helpers
"""    


def get_options():
    parser = OptionParser()

    # Input parameters.
    parser.add_option("-r", "--reads", dest="unpaired_reads_filenames", help="Unpaired FASTQ filenames.")
    parser.add_option("-1", "--first", dest="first_mate_filenames", help="First mate FASTQ filenames.")
    parser.add_option("-2", "--second", dest="second_mate_filenames", help="Second mate FASTQ filenames.")

    # Output parameters.
    parser.add_option("-o", "--output_dir", dest="output_dir", help="Output directory.")

    (options, args) = parser.parse_args()

    return (options,args)


def main():
   
    (options, args) = get_options()

    shell_file = options.output_dir + "/commands.sh"

    ensure_dir(shell_file)
 
    global shell_file_fp
    shell_file_fp = open(shell_file, 'w')
    setup_shell_file()

    # Gather all the filenames.
    reads_filenames = []
    if options.unpaired_reads_filenames:
        reads_filenames.extend(options.unpaired_reads_filenames.split(','))
    if options.first_mate_filenames:
        reads_filenames.extend(options.first_mate_filenames.split(','))
    if options.second_mate_filenames:
        reads_filenames.extend(options.second_mate_filenames.split(','))

    options.reads_filenames = reads_filenames

    # Compress and then decompress the reads.
    compress(options)

    # Assemble the sequences with ALLPATHS-LG.
    assemble(options)



if __name__ == '__main__':
    main()


