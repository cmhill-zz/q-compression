import numpy as np
from itertools import imap
from multiprocessing import Pool
import sys
from array import array

polyfit = np.polynomial.polynomial.polyfit
polyval = np.polyval

def qual_ord(phred_string, offset=33):
    return range(len(phred_string)), map(lambda x: (ord(x)-offset), phred_string)

def map_qualstring(qual_string):
    x, y = qual_ord(qual_string, offset=0)
    coeffs = polyfit(x, y, 2)
    y_out = polyval(coeffs[::-1], x)
    return ''.join([chr(min(max(int(yy),35),74)) for yy in y_out])

def coeffs_for_qualstring(qual_string, degree=5):
    x, y = qual_ord(qual_string, offset=0)
    coeffs = polyfit(x, y, degree)
    return coeffs

def process(lines, pool):
    qual_lines = lines[3::4]
    new_qual_lines = pool.map(map_qualstring, qual_lines)
    lines[3::4] = new_qual_lines
    if lines:
        print '\n'.join(lines)

def process_coeffs_only(lines, pool, float_list):
    qual_lines = lines[3::4]
    qual_coeffs = pool.map(coeffs_for_qualstring, qual_lines)
    float_list.extend(qual_coeffs)


with open(sys.argv[1], 'r') as infile:
    pool = Pool(16)
    
    # use this block for getting the fastq after poly compression - decompression
    # (i.e fastq w/ predicted quality strings, even if we'd never store those)
    lines = []
    for i_counter, line in enumerate(infile):
        if (i_counter % 1000) == 0:
            process(lines, pool)
            lines = []
        lines.append(line.rstrip())
    if len(lines) > 0:
        process(lines, pool)

    # # use this block to write binary file with coefficients
    # float_list = []
    # lines = []
    # for i_counter, line in enumerate(infile):
    #     if (i_counter % 1000) == 0:
    #         process_coeffs_only(lines, pool, float_list)
    #         lines = []
    #     lines.append(line.rstrip())
    # if len(lines) > 0:
    #     process_coeffs_only(lines, pool, float_list)
    # with open('coefficients_binary.bin', 'wb') as output_file:
    #     float_array = array('f', float_list)
    #     float_array.tofile(output_file)

