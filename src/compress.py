#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
from collections import deque
from itertools import islice
#from subprocess import call
import subprocess
from optparse import OptionParser
from tempfile import mkstemp
import glob
import os
import random
import re
import shlex
import shutil
import sys
import tempfile
import time
import resource
import locale
#import file

locale.setlocale(locale.LC_ALL, "C")

FNULL = open('/dev/null', 'w')
base_path = os.path.dirname(sys.argv[0])[:-len('src/')]
dry_run = False

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


def call(call_arr, stdout=sys.stdout, stderr=sys.stderr):

    if not dry_run:
        subprocess.call(call_arr, stdout=stdout, stderr=stderr)


def sort_reads_command(options,reads_filename):
    """
    Sort the incoming FASTQ filename.
    """

    SORT_CMD = "fastq-sort " + reads_filename
    call_arr = SORT_CMD.split()
    output_fp = open(options.output_dir + '/sorted/' + os.path.basename(reads_filename), 'w')
    out_cmd(output_fp.name, FNULL.name, call_arr)
    call(call_arr, stdout=output_fp, stderr=FNULL)
    return output_fp.name


def sort_reads(options):
    """
    Sort the FASTQ reads and update the options accordingly.
    """

    ensure_dir(options.output_dir + '/sorted/')

    if options.unpaired_reads_filenames:
        new_filenames = []
        for reads_filenames in options.unpaired_reads_filenames.split(','):
            new_filenames.append(sort_reads_command(options, reads_filenames))
        options.unpaired_reads_filenames = ','.join(new_filenames)

    if options.first_mate_filenames:
        new_filenames = []
        for reads_filenames in options.first_mate_filenames.split(','):
            new_filenames.append(sort_reads_command(options, reads_filenames))
        options.first_mate_filenames = ','.join(new_filenames)

    if options.second_mate_filenames:
        new_filenames = []
        for reads_filenames in options.second_mate_filenames.split(','):
            new_filenames.append(sort_reads_command(options, reads_filenames))
        options.second_mate_filenames = ','.join(new_filenames)


def compress(options):
    """
    Compress the reads using all methods.
    """

    ensure_dir(options.output_dir + '/original/')
    ensure_dir(options.output_dir + '/goodbad/')
    ensure_dir(options.output_dir + '/maxqual/')
    ensure_dir(options.output_dir + '/minqual/')

    std_err_file = open('compress.log', 'w')

    # Basic command line scripts to run the individual compression schemes.
    GB_COMPRESSION_CMD = "./src/good_bad_coding.py -r [READ] -c 2 -b 0 -i [COMPRESSED_FILE]"
    MAX_VALUE_COMPRESSION_CMD = "./src/good_bad_coding.py -r [READ] -g 40 -b 40 -c 2 -i [COMPRESSED_FILE]"
    MIN_VALUE_COMPRESSION_CMD = "./src/good_bad_coding.py -r [READ] -g 0 -b 0 -c 2 -i [COMPRESSED_FILE]"
    POLY_REGRESSION_CMD = "Rscript src/poly_regression_parallel.R [READ] [OUTPUT] [DEGREE] [COMPRESSED_FILE] [NUM_THREADS]"
    PROFILE_COMPRESSION_CMD = "Rscript src/profile_parallel.R [READ] [OUTPUT] [TRAINING_SIZE] [NUM_PROFILES] [COMPRESSED_FILE] [NUM_THREADS]"

    QUALCOMP_COMPRESS_CMD = "./runCompress.sh -i [READ] -c [CLUSTERS] -r [RATE]"
    QUALCOMP_DECOMPRESS_CMD = "./runDecompress.sh -p [DIR] -c [CLUSTERS] -r [RATE]"

    RQS_COMPRESS_CMD = "./src/run_rqs.sh [READ] [OUTPUT]"

    #qvz -c 3 -r .10 -v test_results4/original/frag_1.fastq.quals tmp/test_1_c3_r.10
    QVZ_COMPRESS_CMD = "[QVZ]/qvz -c [CLUSTERS] -r [RATE] -v [READ] [OUTPUT]"
    QVZ_DECOMPRESS_CMD = "[QVZ]/qvz -x -v [INPUT] [OUTPUT]"

    # Store which compression directories we created.
    options.compressed_dirs = []
    options.compressed_dirs.append('original')
    options.compressed_dirs.append('goodbad')
    options.compressed_dirs.append('maxqual')
    options.compressed_dirs.append('minqual')

    for reads_filename in options.reads_filenames:

        # Copy the original sequences over.
        out_cmd("", std_err_file.name, ["cp", reads_filename, options.output_dir + '/original/' + os.path.basename(reads_filename)])
        shutil.copyfile(reads_filename, options.output_dir + '/original/' + os.path.basename(reads_filename))

        # Good/bad binary compression.
        call_arr = GB_COMPRESSION_CMD.replace('[READ]', reads_filename)\
                .replace('[COMPRESSED_FILE]', options.output_dir + '/goodbad/' + os.path.basename(reads_filename) + '.comp').split()
        output_fp = open(options.output_dir + '/goodbad/' + os.path.basename(reads_filename), 'w')

        out_cmd(options.output_dir + '/goodbad/' + os.path.basename(reads_filename), std_err_file.name, call_arr)
        call(call_arr, stdout=output_fp, stderr=std_err_file)

        # Max/min quality value compression. We can use good_bad.py script to do this.
        call_arr = MAX_VALUE_COMPRESSION_CMD.replace('[READ]', reads_filename)\
                .replace('[COMPRESSED_FILE]', options.output_dir + '/maxqual/' + os.path.basename(reads_filename) + '.comp').split()
        output_fp = open(options.output_dir + '/maxqual/' + os.path.basename(reads_filename), 'w')

        out_cmd(options.output_dir + '/maxqual/' + os.path.basename(reads_filename), std_err_file.name, call_arr)
        call(call_arr, stdout=output_fp, stderr=std_err_file)

        call_arr = MIN_VALUE_COMPRESSION_CMD.replace('[READ]', reads_filename)\
                .replace('[COMPRESSED_FILE]', options.output_dir + '/minqual/' + os.path.basename(reads_filename) + '.comp').split()
        output_fp = open(options.output_dir + '/minqual/' + os.path.basename(reads_filename), 'w')

        out_cmd(options.output_dir + '/minqual/' + os.path.basename(reads_filename), std_err_file.name, call_arr)
        call(call_arr, stdout=output_fp, stderr=std_err_file)

        #continue

        # Polynomial regression.
        if options.poly_degrees:
            for degree in options.poly_degrees.split(','):
                ensure_dir(options.output_dir + '/degree_' + degree + '/')

                if 'degree_' + degree not in options.compressed_dirs:
                    options.compressed_dirs.append('degree_' + degree)

                #continue

                call_arr = POLY_REGRESSION_CMD.replace('[READ]', reads_filename)\
                        .replace('[OUTPUT]', options.output_dir + '/degree_' + degree + '/' + os.path.basename(reads_filename))\
                        .replace('[DEGREE]', degree)\
                        .replace('[COMPRESSED_FILE]', options.output_dir + '/degree_' + degree +'/' + os.path.basename(reads_filename) + '.comp')\
                        .replace('[NUM_THREADS]', options.threads).split()

                out_cmd("", std_err_file.name, call_arr)
                call(call_arr, stderr=std_err_file)

        # Profile compression using k-means.
        if options.profile_sizes:
            for profiles in options.profile_sizes.split(','):
                ensure_dir(options.output_dir + '/profile_' + profiles + '/')

                if 'profile_' + profiles not in options.compressed_dirs:
                    options.compressed_dirs.append('profile_' + profiles)

                #continue

                call_arr = PROFILE_COMPRESSION_CMD.replace('[READ]', reads_filename)\
                        .replace('[OUTPUT]', options.output_dir + '/profile_' + profiles + '/' + os.path.basename(reads_filename))\
                        .replace('[NUM_PROFILES]', profiles)\
                        .replace('[TRAINING_SIZE]', options.training_size)\
                        .replace('[COMPRESSED_FILE]', options.output_dir + '/profile_' + profiles +'/' + os.path.basename(reads_filename) + '.comp')\
                        .replace('[NUM_THREADS]', options.threads).split()

                out_cmd("", std_err_file.name, call_arr)
                call(call_arr, stderr=std_err_file)

        # Compress using QualComp.
        if options.rates:
            for rate in options.rates.split(','):
                #continue
                ensure_dir(options.output_dir + '/qualcomp_r' + rate + '/')

                if 'qualcomp_r' + rate not in options.compressed_dirs:
                    options.compressed_dirs.append('qualcomp_r' + rate)

                #continue

                """
                QUALCOMP_COMPRESS_CMD = "$QUALCOMP/runCompressMod.sh -i [READ] -c [CLUSTERS] -r [RATE]"
                QUALCOMP_DECOMPRESS_CMD = "$QUALCOMP/runDecompress.sh -p [DIR] -c [CLUSTERS] -r [RATE]"
                """

                reads_abs_path = os.path.abspath(reads_filename)
                prev_dir = os.getcwd()
                os.chdir(os.environ["QUALCOMP"])

                call_arr = QUALCOMP_COMPRESS_CMD.replace('[READ]', reads_abs_path)\
                        .replace('[CLUSTERS]', options.clusters)\
                        .replace('[RATE]', rate).split()

                out_cmd(std_err_file.name, std_err_file.name, call_arr)
                call(call_arr, stdout=std_err_file, stderr=std_err_file)

                # Also decompress using QualComp special function.
                qualcomp_prefix = reads_abs_path.split('.')[0]
                call_arr = QUALCOMP_DECOMPRESS_CMD.replace('[DIR]', qualcomp_prefix)\
                        .replace('[CLUSTERS]', options.clusters)\
                        .replace('[RATE]', rate).split()

                out_cmd(std_err_file.name, std_err_file.name, call_arr)
                call(call_arr, stdout=std_err_file, stderr=std_err_file)

                os.chdir(prev_dir)

                # QualComp writes the files into the original directory,
                # so move the fastq files into the QualComp directory.
                mv_cmd = "mv " + qualcomp_prefix + "_" + options.clusters + "_" + rate + ".fastq " + options.output_dir + '/qualcomp_r' + rate + '/' + os.path.basename(reads_filename)
                call_arr = mv_cmd.split()
                out_cmd("", std_err_file.name, call_arr)
                call(call_arr, stderr=std_err_file)

                filename_list = glob.glob(qualcomp_prefix + "_" + options.clusters + "_*")
                mv_cmd = "mv " + ' '.join(filename_list) + ' ' + options.output_dir + '/qualcomp_r' + rate + '/'
                call_arr = mv_cmd.split()
                out_cmd("", std_err_file.name, call_arr)
                call(call_arr, stderr=std_err_file)

                # Concatenate all the binary files to create a single 'compressed' file.
                filename_list = glob.glob(options.output_dir + '/qualcomp_r' + rate + '/' + os.path.basename(reads_filename).split(".")[0] +  "*bin")
                cat_cmd = "cat " + ' '.join(filename_list)
                call_arr = cat_cmd.split()
                bin_file = open(options.output_dir + '/qualcomp_r' + rate + '/' + os.path.basename(reads_filename) + '.comp', 'w')
                out_cmd(bin_file.name, std_err_file.name, call_arr)
                call(call_arr, stdout=bin_file, stderr=std_err_file)

        # Compress using QualComp.
        if options.qvz_rates:
            for rate in options.qvz_rates.split(','):
                ensure_dir(options.output_dir + '/qvz_r' + rate + '/')

                if 'qvz_r' + rate not in options.compressed_dirs:
                    options.compressed_dirs.append('qvz_r' + rate)

                """
                QVZ_COMPRESS_CMD = "[QVZ]/qvz -c [CLUSTERS] -r [RATE] -v [READ] [OUTPUT]"
                QVZ_DECOMPRESS_CMD = "[QVZ]/qvz -x -v [INPUT] [OUTPUT]"
                """

                # We first need to get the quality values only. TODO(cmhill): Very hacky.
                with open(reads_filename) as fin, open(options.output_dir + '/qvz_r' + rate + '/orig_' +  os.path.basename(reads_filename) + '.quals', 'w') as fout:
                    fout.writelines(islice(fin, 3, None, 4))

                call_arr = QVZ_COMPRESS_CMD.replace('[READ]', options.output_dir + '/qvz_r' + rate + '/orig_' +  os.path.basename(reads_filename) + '.quals')\
                        .replace('[QVZ]', os.environ["QVZ"])\
                        .replace('[CLUSTERS]', options.clusters)\
                        .replace('[RATE]', rate)\
                        .replace('[OUTPUT]', options.output_dir + '/qvz_r' + rate + '/' + os.path.basename(reads_filename) + '.comp').split()

                out_cmd(std_err_file.name, std_err_file.name, call_arr)
                call(call_arr, stdout=std_err_file, stderr=std_err_file)

                # Also decompress using QVZ special function.
                #qualcomp_prefix = reads_abs_path.split('.')[0]
                call_arr = QVZ_DECOMPRESS_CMD.replace('[INPUT]', options.output_dir + '/qvz_r' + rate + '/' + os.path.basename(reads_filename) + '.comp')\
                        .replace('[QVZ]', os.environ["QVZ"])\
                        .replace('[OUTPUT]', options.output_dir + '/qvz_r' + rate + '/' + os.path.basename(reads_filename) + '.quals').split()

                out_cmd(std_err_file.name, std_err_file.name, call_arr)
                call(call_arr, stdout=std_err_file, stderr=std_err_file)

                # Rebuild the FASTQ file from the unconstructed quals.
                call_arr_str = "python ./src/merge_qual_into_fastq.py " + reads_filename + " " + options.output_dir + '/qvz_r' + rate + '/' + os.path.basename(reads_filename) + '.quals'
                call_arr = call_arr_str.split()
                decompressed_file = open(options.output_dir + '/qvz_r' + rate + '/' + os.path.basename(reads_filename), 'w')
                out_cmd(decompressed_file.name, std_err_file.name, call_arr)
                call(call_arr, stdout=decompressed_file, stderr=std_err_file)

        # Compress with RQS.
        ensure_dir(options.output_dir + '/rqs/')

        if 'rqs' not in options.compressed_dirs:
            options.compressed_dirs.append('rqs')

        call_arr = RQS_COMPRESS_CMD.replace('[READ]', reads_filename)\
                .replace('[OUTPUT]', options.output_dir + '/rqs/').split()

        out_cmd("", std_err_file.name, call_arr)
        call(call_arr, stderr=std_err_file)

    # # Profile compression using k-means.
    # for profiles in options.profile_sizes.split(','):
    #     ensure_dir(options.output_dir + '/profile_' + profiles + '/')
    #
    #     if 'profile_' + profiles not in options.compressed_dirs:
    #         options.compressed_dirs.append('profile_' + profiles)
    #
    #     #continue
    #
    #     call_arr = PROFILE_COMPRESSION_CMD.replace('[READ]', reads_filename)\
    #             .replace('[OUTPUT]', options.output_dir + '/profile_' + profiles + '/' + os.path.basename(reads_filename))\
    #             .replace('[NUM_PROFILES]', profiles)\
    #             .replace('[TRAINING_SIZE]', options.training_size)\
    #             .replace('[COMPRESSED_FILE]', options.output_dir + '/profile_' + profiles +'/' + os.path.basename(reads_filename) + '.comp')\
    #             .replace('[NUM_THREADS]', options.threads).split()
    #
    #     out_cmd("", std_err_file.name, call_arr)
    #     call(call_arr, stderr=std_err_file)

    # After we compress/decompress everything, write out the quality values to a separate file and then run bzip on them.
    for compression_method in options.compressed_dirs:
        for reads_filename in options.reads_filenames:

            decompressed_file = options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename)

            # Check if we already made a quals file, like in the case of QVZ.
            if not os.path.isfile(decompressed_file + '.quals'):
                with open(decompressed_file) as fin, open(decompressed_file + '.quals', 'w') as fout:
                    fout.writelines(islice(fin, 3, None, 4))

            # Even though we do it in python, output the awk command in case someone runs it independently.
            cmd = 'awk \'{if (NR % 4 == 0) print $0}\' ' + options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename)
            out_cmd(decompressed_file + '.quals', std_err_file.name, 'awk \'{if (NR % 4 == 0) print $0}\''.split())

            # Bzip2 the quality values.
            cmd = "bzip2 -k " + options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.quals'
            out_cmd("", std_err_file.name, cmd.split())
            call(cmd.split(), stderr=std_err_file)

            # Bzip2 the compressed quality values.
            if os.path.isfile(options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.comp'):
                cmd = "bzip2 -k " + options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.comp'
                out_cmd("", std_err_file.name, cmd.split())
                call(cmd.split(), stderr=std_err_file)


    # Calculate the information lost from compression.
    calc_mean_squared_error(options)

    std_err_file.close()

    pass


def calc_mean_squared_error(options):
    """
    Calculate mean squared error between the original and decompressed reads.
    """

    std_err_file = open('mse.log', 'w')

    for compression_method in options.compressed_dirs:

        for reads_filename in options.reads_filenames:
            MSE_CMD = "python src/evaluate_loss.py " + options.output_dir + '/original/' + os.path.basename(reads_filename) + '.quals' + '\t' + \
                    options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.quals'
            mse_std_out = open(options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.mse', 'w')
            out_cmd(mse_std_out.name, std_err_file.name, MSE_CMD.split())
            call(MSE_CMD.split(), stdout=mse_std_out, stderr=std_err_file)


def quality_preprocessing(options):
    """
    Examine the effects of lossy compression on quality preprocessing tools.
    """
    std_err_file = open('preprocessing.log', 'w')

    SICKLE_CMD = "sickle se -f [READ] -t sanger -o [OUTPUT]"

    for compression_method in options.compressed_dirs:

        for reads_filename in options.reads_filenames:
            output_filename = options.output_dir + '/preprocessing/' + compression_method + '/' + os.path.basename(reads_filename)
            ensure_dir(output_filename)

            stats_file = open(options.output_dir + '/preprocessing/' + compression_method + '/' + os.path.basename(reads_filename) + '.stats', 'w')

            call_arr = SICKLE_CMD.replace('[READ]', options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename)).replace('[OUTPUT]', output_filename).split()
            out_cmd("", stats_file.name, call_arr)
            call(call_arr, stdout=stats_file)

            # Process the stats file to get the relevant information.
            """Output stats file is in the format of:

            SE input file: tmp/goodbad/frag_small_2.fastq

            Total FastQ records: 200
            FastQ records kept: 138
            FastQ records discarded: 62
            """

            # Parse the above to get how many records were kept.
            for line in open(stats_file.name, 'r').readlines():
                if line.startswith('FastQ records kept:'):
                    records_kept = line.strip().split()[3]
                    open(stats_file.name + '.records_kept', 'w').write(records_kept + '\n')

            # Find out many bases were kept and write out the header files to a separate file.
            line_number = 1
            bases = 0
            headers = []
            for line in open(output_filename, 'r'):
                if (line_number % 4) == 1:
                    headers.append(line.strip())
                if (line_number % 4) == 2:
                    bases += len(line.strip())
                line_number += 1

            open(stats_file.name + '.bases', 'w').write(str(bases) + '\n')

            headers.sort(cmp=locale.strcoll)
            open(options.output_dir + '/preprocessing/' + compression_method + '/' + os.path.basename(reads_filename) + '.headers', 'w').write('\n'.join(headers) + '\n')



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
        call(call_arr, stderr=std_err_file)

        # For AllpathsLG to run successfully, need to have a PLOIDY file present.
        PLOIDY_CMD = "echo 1"
        call_arr = PLOIDY_CMD.split()
        ploidy_file = open(options.output_dir + '/assemble/' + compression_method + '/ploidy', 'w')
        out_cmd(ploidy_file.name, std_err_file.name, call_arr)
        call(call_arr, stdout=ploidy_file, stderr=std_err_file)

        # Run AllpathsLG
        ALLPATHS_CMD = "RunAllPathsLG PRE=" + os.path.abspath(options.output_dir) + '/assemble/' + compression_method + " DATA_SUBDIR=. RUN=allpaths SUBDIR=run THREADS=" + options.threads + " OVERWRITE=True REFERENCE_NAME=."
        call_arr = ALLPATHS_CMD.split()
        out_cmd("", std_err_file.name, call_arr)
        call(call_arr, stderr=std_err_file)

        # Calculate the assembly likelihood.
        # /cbcb/project-scratch/cmhill/metalap/calc_prob.py -1 original/frag_1.fastq -2 original/frag_2.fastq -q -a decomp_0/allpaths/ASSEMBLIES/run/final.assembly.fasta  -I 0 -X 500  -m 180 -t 18  -p 32
        reads_filenames = ','.join(options.reads_filenames)
        assembly_filename = os.path.abspath(options.output_dir + '/assemble/' + compression_method + '/allpaths/ASSEMBLIES/run/final.assembly.fasta')
        CALC_PROB_CMD = os.environ["CALCPROB"] + "/calc_prob.py -q -i " + reads_filenames + ' -a ' + assembly_filename + ' -p ' + options.threads

        # tmp_rhodo/assemble/original/allpaths/ASSEMBLIES/run/final.assembly.fasta
        probs_file = open(options.output_dir + '/assemble/' + compression_method + '/output.probs', 'w')

        call_arr = CALC_PROB_CMD.split()
        out_cmd(probs_file.name, std_err_file.name, call_arr)
        call(call_arr, stdout=probs_file, stderr=std_err_file)

        # Sum the read probabilities.
        SUM_PROB_CMD = os.environ["CALCPROB"] + '/sum_prob.py -t 1e-80 -i ' + probs_file.name
        sum_file = open(options.output_dir + '/assemble/' + compression_method + '/output.sum', 'w')
        call_arr = SUM_PROB_CMD.split()
        out_cmd(sum_file.name, std_err_file.name, call_arr)
        call(call_arr, stdout=sum_file, stderr=std_err_file)

        # Run the getCorrectnessStatistics.sh
        #  sh getCorrectnessStats.sh [reference] [assembly] [assembly]
        #reads_abs_path = os.path.abspath(reads_filename)
        output_abs_dir = os.path.abspath(options.output_dir + '/assemble/' + compression_method + '/')
        reference_abs_path = os.path.abspath(options.reference_fasta)
        prev_dir = os.getcwd()
        os.chdir(os.environ["MUMMER"])
        GET_CORRECTNESS_STATS_CMD = "sh ./getCorrectnessStats.sh " + reference_abs_path + ' ' + assembly_filename + ' ' + assembly_filename
        correctness_stats_file = open(output_abs_dir + '/assembly.correctness', 'w')
        call_arr = GET_CORRECTNESS_STATS_CMD.split()
        out_cmd(correctness_stats_file.name, std_err_file.name, call_arr)
        call(call_arr, stdout=correctness_stats_file, stderr=std_err_file)

        os.chdir(prev_dir)


def align_reads(options):
    """
    Evaluate how all decompressed reads align with Bowtie2.
    """

    std_err_file = open('alignment.log', 'a')

    # Construct the Bowtie2 index.
    ensure_dir(options.output_dir + "/align/")
    BOWTIE2_INDEX_CMD = "bowtie2-build " + options.reference_fasta +  " " + options.output_dir + "/align/reference"
    call_arr = BOWTIE2_INDEX_CMD.split()
    out_cmd("", "", call_arr)
    call(call_arr, stdout=FNULL, stderr=std_err_file)

    # Align the reads.
    BOWTIE2_CMD = "bowtie2 -x " + options.output_dir + "/align/reference -p " + options.threads + " --reorder -U [READ] --al [ALIGNED] -S [SAM]"

    # Have to do max/min alignment

    for compression_method in options.compressed_dirs:
        for reads_filename in options.reads_filenames:

            alignment_filename = options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.alignment_summary'
            ensure_dir(alignment_filename)
            alignment_file = open(alignment_filename, 'w')

            call_arr = BOWTIE2_CMD.replace('[READ]', options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename))\
                    .replace('[ALIGNED]', options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.aligned')\
                    .replace('[SAM]', options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.sam').split()
            out_cmd(FNULL.name, alignment_filename, call_arr)
            call(call_arr, stdout=FNULL, stderr=alignment_file)

            # Print out the headers of the aligned sequences.
            with open(options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.aligned', 'r') as aligned_file, open(options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.headers', 'w') as headers:
                headers.writelines(islice(aligned_file, 0, None, 4))

def call_snps(options):
    """
    Call SNPs using decompressed reads.
    """

    std_err_file = open('snp.log', 'a')

    ensure_dir(options.output_dir + "/snp/")
    SAMTOOLS_SORT_CMD = "samtools sort -T ./ -o [SORTED_SAM] -O sam [SAM]"
    PILEUP_CMD = "samtools mpileup -ugf [REFERENCE] [SAM]"
    SNP_CMD = "bcftools call -vV indels -mO v -o [OUTPUT_VCF] [PILEUP]"

    for compression_method in options.compressed_dirs:
        for reads_filename in options.reads_filenames:

            pileup_filename = options.output_dir + '/snp/' + compression_method + '/' + os.path.basename(reads_filename) + '.pileup'
            ensure_dir(pileup_filename)

            call_arr = SAMTOOLS_SORT_CMD.replace('[SAM]', options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.sam')\
                    .replace('[SORTED_SAM]', options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.sorted.sam').split()
            out_cmd(std_err_file.name, std_err_file.name, call_arr)
            call(call_arr, stdout=std_err_file, stderr=std_err_file)

            pileup_file = open(pileup_filename, 'w')
            call_arr = PILEUP_CMD.replace('[REFERENCE]', options.reference_fasta)\
                    .replace('[SAM]', options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.sorted.sam').split()
            out_cmd(pileup_file.name, std_err_file.name, call_arr)
            call(call_arr, stdout=pileup_file, stderr=std_err_file)

            call_arr = SNP_CMD.replace('[OUTPUT_VCF]', options.output_dir + '/snp/' + compression_method + '/' + os.path.basename(reads_filename) + '.vcf')\
                    .replace('[PILEUP]', options.output_dir + '/snp/' + compression_method + '/' + os.path.basename(reads_filename) + '.pileup').split()
            out_cmd(std_err_file.name, std_err_file.name, call_arr)
            call(call_arr, stdout=std_err_file, stderr=std_err_file)

            # Print out the locations of the SNP.
            snps_filename = options.output_dir + '/snp/' + compression_method + '/' + os.path.basename(reads_filename) + '.snps'
            with open(snps_filename, 'w') as snps_file, open(options.output_dir + '/snp/' + compression_method + '/' + os.path.basename(reads_filename) + '.vcf', 'r') as vcf_file:
                line = vcf_file.readline()
                while line.startswith('#'):
                    line = vcf_file.readline()

                snps = []
                while line:
                    tuple = line.split('\t')
                    if len(tuple)>1:
                        snps.append(str(tuple[0] + '\t' + tuple[1] + '\n'))
                        #snps_file.write(str(tuple[0] + '\t' + tuple[1] + '\n'))

                    line = vcf_file.readline()

                snps.sort()
                for snp in snps:
                    snps_file.write(snp)


def post_process_results(options):
    """
    Consolidate the results.
    """

    process_compression_stats(options)

    if options.preprocessing:
        process_preprocessing_stats(options)

    if options.assemble:
        process_assembly_stats(options)

    if options.alignment or options.snp:
        process_alignment_stats(options)

        if options.snp:
            process_snp_stats(options)



def process_compression_stats(options):
    """
    Ouput compression results for each FASTQ file:
            compression_method original_size compressed_size bzip2_size
    """

    ensure_dir(options.output_dir + "/results/")

    for reads_filename in options.reads_filenames:

        # Store the amount of bases we see.
        bases = 0
        for line in open(options.output_dir + '/original/' + os.path.basename(reads_filename) + '.quals', 'r'):
            bases += len(line.strip())


        results_file = open(options.output_dir + "/results/" + os.path.basename(reads_filename) + '.compression', 'w')
        results_file.write("compression\tbases\torig_size\torig_bzip2\tcomp_size\tcomp_bzip2\tmse\tbits_bp\n")
        for compression_method in options.compressed_dirs:

            filename = options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename)

            results = compression_method + '\t'

            results += str(bases) + '\t'

            # Get the original size.
            results += str(os.path.getsize(filename + '.quals')) + '\t'

            # Get the bzip2 size.
            results += str(os.path.getsize(filename + '.quals.bz2')) + '\t'

            compressed_size = os.path.getsize(filename + '.quals.bz2')

            # Get the compressed size.
            if os.path.isfile(filename + '.comp'):
                results += str(os.path.getsize(filename + '.comp')) + '\t'
                results += str(os.path.getsize(filename + '.comp.bz2')) + '\t'
                if os.path.getsize(filename + '.comp.bz2') < compressed_size:
                    compressed_size = os.path.getsize(filename + '.comp.bz2')
            else:
                results += "NA\tNA\t"

            # Get the MSE.
            results += str(grab_value_from_file(filename + '.mse')) + '\t'

            # Print the bits/bp.
            results += str((compressed_size * 8) / float(bases)) + '\n'

            results_file.write(results)

        results_file.close()


def process_preprocessing_stats(options):
    """
    Output the preprocessing results.
    """

    for reads_filename in options.reads_filenames:
        sequence_count = num_lines(reads_filename) / 4

        results_file = open(options.output_dir + "/results/" + os.path.basename(reads_filename) + '.preprocessing', 'w')
        results_file.write("compression\tuniq_orig\tuniq_comp\tcommon_kept\tcommon_discard\tbases_kept\n")
        for compression_method in options.compressed_dirs:

            results = compression_method + '\t'

            # Get the number of sequences kept unique to original.
            unique_to_orig = run_comm_and_return_line_count(options, "23", options.output_dir + '/preprocessing/original/'\
                    + os.path.basename(reads_filename) + '.headers',\
                    options.output_dir + '/preprocessing/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.headers')
            results += unique_to_orig + '\t'

            # Get the number of sequences kept unique to the compression method.
            unique_to_compress = run_comm_and_return_line_count(options, "13", options.output_dir + '/preprocessing/original/'\
                    + os.path.basename(reads_filename) + '.headers',\
                    options.output_dir + '/preprocessing/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.headers')
            results += unique_to_compress + '\t'

            # Get the number of sequences kept by both methods.
            common_count = run_comm_and_return_line_count(options, "12", options.output_dir + '/preprocessing/original/'\
                    + os.path.basename(reads_filename) + '.headers',\
                    options.output_dir + '/preprocessing/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.headers')
            results += common_count + '\t'

            # Get the number of sequences filtered by both methods.
            results += str(sequence_count - int(unique_to_orig) - int(unique_to_compress) - int(common_count)) + '\t'

            # Get the bases kept in the file.
            results += str(grab_value_from_file(options.output_dir + '/preprocessing/' + compression_method + '/' + os.path.basename(reads_filename) + '.stats.bases'))

            results += '\n'
            results_file.write(results)

        results_file.close()


def process_assembly_stats(options):
    """
    Output the assembly stats.
    """

    results_file = open(options.output_dir + "/results/assembly", 'w')
    results_file.write("compression\tref_bases\tasm_bases\tN50\tmissing_ref_bases\tavg_idy\tsnps\tindels_gt5bp\tinversions\treloc\ttransloc\tcorr_N50\tLAP\n")

    for compression_method in options.compressed_dirs:
        results = compression_method + '\t'

        results += "\t".join(parse_correctness_stats(options.output_dir + '/assemble/' + compression_method + '/assembly.correctness'))

        results += '\t' + str(grab_value_from_file(options.output_dir + '/assemble/' + compression_method + '/output.sum'))

        results += '\n'
        results_file.write(results)

    results_file.close()


def parse_correctness_stats(filename):
    """
    Parse the results returned from get_Correctness.sh
    """

    results = []

    stats_file = open(filename, 'r')

    line = stats_file.readline()
    while not line.startswith('Reference:'):
        line = stats_file.readline()

    # Add Reference bases.
    results.append(line.strip().split()[1])

    # Add the genome bases.
    line = stats_file.readline()
    results.append(line.strip().split()[1])

    line = stats_file.readline()
    while not line.startswith('N50:'):
        line = stats_file.readline()

    # Add N50.
    results.append(line.split()[1])

    line = stats_file.readline()
    while not line.startswith('Missing Reference'):
        line = stats_file.readline()

    # Missing reference bases.
    results.append(line.strip().split(' ')[3].split('(')[0])

    line = stats_file.readline()
    while not line.startswith('Avg Idy:'):
        line = stats_file.readline()

    # Add Avg Idy.
    results.append(line.strip().split(' ')[2])

    # Add SNPs.
    line = stats_file.readline()
    results.append(line.strip().split(' ')[1])

    # Add Indels > 5bp.
    line = stats_file.readline()
    line = stats_file.readline()
    results.append(line.strip().split(' ')[3])

    # Add Inversions.
    line = stats_file.readline()
    results.append(line.strip().split(' ')[1])

    # Add Relocations.
    line = stats_file.readline()
    results.append(line.strip().split(' ')[1])

    # Add translocations.
    line = stats_file.readline()
    results.append(line.strip().split(' ')[1])

    line = stats_file.readline()
    while not line.startswith('N50:'):
        line = stats_file.readline()

    # Add Corrected N50.
    results.append(line.strip().split()[1])

    stats_file.close()

    return results


def process_alignment_stats(options):
    """
    Process the bowtie2 alignment results.
    """

    # options.output_dir + '/' + compression_method + '/' + os.path.basename(reads_filename) + '.aligned'
    for reads_filename in options.reads_filenames:
        sequence_count = num_lines(reads_filename) / 4

        results_file = open(options.output_dir + "/results/" + os.path.basename(reads_filename) + '.alignment', 'w')
        results_file.write("compression\tmapped\tshared\tuniq_orig\tuniq_comp\tunmapped_shared\n")
        for compression_method in options.compressed_dirs:

            results = compression_method + '\t'

            # Get aligned reads.
            results += str(num_lines(options.output_dir + '/align/' + compression_method + '/'\
                    + os.path.basename(reads_filename) + '.headers')) + '\t'

            # Get the number of sequences kept by both methods.
            common_count = run_comm_and_return_line_count(options, "12", options.output_dir + '/align/original/'\
                    + os.path.basename(reads_filename) + '.headers',\
                    options.output_dir + '/align/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.headers')
            results += common_count + '\t'

            # Get the number of sequences kept unique to original.
            unique_to_orig = run_comm_and_return_line_count(options, "23", options.output_dir + '/align/original/'\
                    + os.path.basename(reads_filename) + '.headers',\
                    options.output_dir + '/align/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.headers')
            results += unique_to_orig + '\t'

            # Get the number of sequences kept unique to the compression method.
            unique_to_compress = run_comm_and_return_line_count(options, "13", options.output_dir + '/align/original/'\
                    + os.path.basename(reads_filename) + '.headers',\
                    options.output_dir + '/align/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.headers')
            results += unique_to_compress + '\t'

            # Get the number of sequences filtered by both methods.
            results += str(sequence_count - int(unique_to_orig) - int(unique_to_compress) - int(common_count))

            # Get the bases kept in the file.
            #results += str(grab_value_from_file(options.output_dir + '/align/' + compression_method + '/' + os.path.basename(reads_filename) + '.stats.bases'))

            results += '\n'
            results_file.write(results)

        results_file.close()


def process_snp_stats(options):
    """
    Process the SNP results.
    """

    for reads_filename in options.reads_filenames:
        sequence_count = num_lines(reads_filename) / 4

        results_file = open(options.output_dir + "/results/" + os.path.basename(reads_filename) + '.snps', 'w')
        results_file.write("compression\tshared\tuniq_orig\tuniq_comp\n")
        for compression_method in options.compressed_dirs:
            results = compression_method + '\t'

            # Get the number of SNPs kept by both methods.
            common_count = run_comm_and_return_line_count(options, "12", options.output_dir + '/snp/original/'\
                    + os.path.basename(reads_filename) + '.snps',\
                    options.output_dir + '/snp/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.snps')
            results += common_count + '\t'

            # Get the number of SNPs unique to original.
            unique_to_orig = run_comm_and_return_line_count(options, "23", options.output_dir + '/snp/original/'\
                    + os.path.basename(reads_filename) + '.snps',\
                    options.output_dir + '/snp/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.snps')
            results += unique_to_orig + '\t'

            # Get the number of SNPs unique to the compression method.
            unique_to_compress = run_comm_and_return_line_count(options, "13", options.output_dir + '/snp/original/'\
                    + os.path.basename(reads_filename) + '.snps',\
                    options.output_dir + '/snp/' + compression_method + '/'\
                     + os.path.basename(reads_filename) + '.snps')
            results += unique_to_compress

            results += '\n'
            results_file.write(results)

        results_file.close()


def run_comm_and_return_line_count(options, suppress, file1, file2):
    """
    Run: comm -1 [suppress] [file1] [file2] | wc -l
    """

    std_err_file = open('preprocessing.log', 'a')

    cmd = "comm -" + suppress + " " + file1 + " " + file2
    tmp_file = tempfile.NamedTemporaryFile('w', dir=options.output_dir)

    call_arr = cmd.split()
    out_cmd(tmp_file.name, "", call_arr)
    call(call_arr, stdout=tmp_file, stderr=std_err_file)


    cmd = "wc -l " + tmp_file.name
    tmp_wc_file = tempfile.NamedTemporaryFile('w', dir=options.output_dir)
    call_arr = cmd.split()
    out_cmd(tmp_wc_file.name, "", call_arr)
    call(call_arr, stdout=tmp_wc_file, stderr=std_err_file)

    count = grab_value_from_file(tmp_wc_file.name)

    tmp_file.close()
    tmp_wc_file.close()

    return count


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


def grab_value_from_file(filename):
    """
    Return the value from the first line of a file.
    """
    return open(filename, 'r').readline().strip().split()[0]

def num_lines(filename):
    """
    http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    """
    return sum(1 for line in open(filename, 'r'))

"""
I/O Helpers
"""


def get_options():
    parser = OptionParser()

    # Input parameters.
    parser.add_option("-r", "--reads", dest="unpaired_reads_filenames", help="Unpaired FASTQ filenames.")
    parser.add_option("-1", "--first", dest="first_mate_filenames", help="First mate FASTQ filenames.")
    parser.add_option("-2", "--second", dest="second_mate_filenames", help="Second mate FASTQ filenames.")
    parser.add_option("-f", "--reference", dest="reference_fasta", help="Reference FASTA filename.")

    # Output parameters.
    parser.add_option("-o", "--output_dir", dest="output_dir", help="Output directory.")

    # Pipeline options.
    parser.add_option("-a", "--assemble", dest="assemble", help="Run assembly evaluation.", action='store_true')
    parser.add_option("-p", "--preprocessing", dest="preprocessing", help="Run preprocessing tools evaluation.", action='store_true')
    parser.add_option("-b", "--alignment", dest="alignment", help="Run alignment evaluation (using Bowtie2).", action='store_true')
    parser.add_option("-n", "--snp", dest="snp", help="Run SNP analysis (will include alignment option if not already selected).", action='store_true')
    parser.add_option("-s", "--sort", dest="sort_reads", help="Sort FASTQ reads before the pipeline begins (requires fastq-sort).", action='store_true')

    # Polynomial regression specific options.
    parser.add_option("--poly-degrees", dest="poly_degrees", help="Comma-separated list of polynomial degrees to use for regression.")

    # Profile-specific compression options.
    parser.add_option("--training-size", dest="training_size", help="Training size used for clustering.", default = "10000")
    parser.add_option("--profile-sizes", dest="profile_sizes", help="Comma-separated list of number of profiles to use.", default="256")

    # QualComp-specific compression options.
    parser.add_option("--rates", dest="rates", help="QualComp parameter for setting the bits/reads.", default=None)#"30")
    parser.add_option("--clusters", dest="clusters", help="QualComp parameter for setting number of clusters.", default="3")

    # RQS-specific compression options.
    #parser.add_option("--rqs", dest="rates", help="QualComp parameter for setting the  bits/reads.", default="30")
    #parser.add_option("--clusters", dest="clusters", help="QualComp parameter for setting number of clusters.", default="3")

    # QVZ-specific compression options.
    parser.add_option("--qvz-rates", dest="qvz_rates", help="QVZ parameter for setting the bits/reads.", default=None)#".30")
    parser.add_option("--qvz-clusters", dest="qvz_clusters", help="QVZ parameter for setting number of clusters.", default="3")

    # Max, min quality value compression options.
    parser.add_option("--max-qv", dest="max_quality", help="Use this value for max quality value compression.", default="40")
    parser.add_option("--min-qv", dest="min_quality", help="Use this value for min quality value compression.", default="10")

    # Additional options.
    parser.add_option("-t", "--threads", dest="threads", help="Number of threads (default 32).", default="32")
    parser.add_option("-d", "--dry-run", dest="dry_run", help="Don't run any commands.", action='store_true')

    (options, args) = parser.parse_args()

    return (options,args)


def main():

    (options, args) = get_options()

    shell_file = options.output_dir + "/commands.sh"

    ensure_dir(shell_file)

    global shell_file_fp
    shell_file_fp = open(shell_file, 'w')
    setup_shell_file()

    global dry_run
    dry_run = options.dry_run

    # Sort the fastq files first if necessary.
    if options.sort_reads:
        sort_reads(options)

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

    # Carry out preprocessing evaluation with SICKLE.
    if options.preprocessing:
        quality_preprocessing(options)

    # Align the reads using Bowtie2.
    if options.alignment or options.snp:
        align_reads(options)

        if options.snp:
            call_snps(options)

    # Output the compression results.
    post_process_results(options)


if __name__ == '__main__':
    main()
