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


    # After we compress/decompress everything, write out the quality values to a separate file and then run bzip on them.
    for compression_method in options.compressed_dirs:
        for reads_filename in options.reads_filenames:
            from itertools import islice

            decompressed_file = options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename)

            with open(decompressed_file) as fin, open(decompressed_file + '.quals', 'w') as fout:
                fout.writelines(islice(fin, 3, None, 4))
            
            # Even though we do it in python, output the awk command in case someone runs it independently.
            cmd = 'awk \'{if (NR % 4 == 0) print $0}\' ' + options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename)
            out_cmd(decompressed_file + '.quals', std_err_file.name, 'awk \'{if (NR % 4 == 0) print $0}\''.split())
            
            # Bzip2 the quality values.
            cmd = "bzip2 -k " + options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.quals'
            out_cmd("", std_err_file.name, cmd.split())
            call(cmd.split(), stderr=std_err_file)


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

    for compression_method in options.compressed_dirs:
        for reads_filename in options.reads_filenames:
            pass
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

    std_err_file = open('assemble.log', 'a')

    # The first thing step is to create the in_groups.csv and in_libs.csv.

    IN_GROUPS_CSV = """group_name,\tlibrary_name,\tfile_name
frag,\tIllumina_01,\t[FULL_PATH]/[COMPRESSION]/frag_*.fastq
shortjump,\tIllumina_02,\t[FULL_PATH]/[COMPRESSION]/shortjump_*.fastq"""

    IN_LIBS_CSV = """library_name,\tproject_name,\torganism_name,\ttype,\tpaired,\tfrag_size,\tfrag_stddev,\tinsert_size,\tinsert_stddev,\tread_orientation,\tgenomic_start,\tgenomic_end
Illumina_01,\tassembly,\tunknown,\tfragment,\t1,\t180,\t10,\t,\t,\tinward,\t,\t
Illumina_02,\tassembly,\tunknown,\tjumping,\t1,\t,\t,\t3000,\t500,\toutward,\t,\t"""

    #print(IN_GROUPS_CSV.replace('[FULL_PATH]', os.path.abspath(options.output_dir)).replace('[COMPRESSION]', 'goodbad'))

    for compression_method in options.compressed_dirs:
        ensure_dir(options.output_dir + '/assemble/' + compression_method + '/')

        open(options.output_dir + '/assemble/' + compression_method + '/in_groups.csv', 'w').write(IN_GROUPS_CSV.replace('[FULL_PATH]', \
                os.path.abspath(options.output_dir)).replace('[COMPRESSION]', compression_method))
        open(options.output_dir + '/assemble/' + compression_method + '/in_libs.csv', 'w').write(IN_LIBS_CSV)

        # Prepare the input for AllpathsLG.
        PREPARE_CMD = 'PrepareAllPathsInputs.pl DATA_DIR=' + os.path.abspath(options.output_dir) + '/assemble/' + compression_method + \
                ' IN_GROUPS_CSV=' + os.path.abspath(options.output_dir) + '/assemble/' + compression_method + '/in_groups.csv' + \
                ' IN_LIBS_CSV=' + os.path.abspath(options.output_dir) + '/assemble/' + compression_method + '/in_libs.csv'

        call_arr = PREPARE_CMD.split()
        out_cmd("", std_err_file.name, call_arr)

        # Run AllpathsLG
        ALLPATHS_CMD = "RunAllPathsLG PRE=" + os.path.abspath(options.output_dir) + '/assemble/' + compression_method + " DATA_SUBDIR=. RUN=allpaths SUBDIR=run THREADS=32 OVERWRITE=True"
        # RunAllPathsLG PRE=. REFERENCE_NAME=. DATA_SUBDIR=. RUN=allpaths SUBDIR=run THREADS=32 OVERWRITE=True
        # RunAllPathsLG PRE=/assemblies DATA=datadir RUN=allpaths SUBDIR=run THREADS=32 OVERWRITE=True

        call_arr = ALLPATHS_CMD.split()
        #output_fp = open(options.output_dir + '/goodbad/' + os.path.basename(reads_filename), 'w')

        out_cmd("", std_err_file.name, call_arr)
        call(call_arr, stderr=std_err_file)


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

    # Pipeline options.
    parser.add_option("-a", "--assemble", dest="assemble", help="Run assembly evaluation", action='store_true')

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
    if options.assemble:
        assemble(options)



if __name__ == '__main__':
    main()


