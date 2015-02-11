#! /usr/bin/env python

"""

"""

from optparse import OptionParser
import math, os, sys

cutoff = 0
good_value = 40
bad_value = 10
offset_value = 33
bin_size = 10

def main():
    
    (options, args) = get_options()
    
    # Open file and store quality values.
    reads_file = open(options.reads_filename, 'r')

    quality_values = []

    """
    @SRR081522.906/1
    NCTGATCCTCGGGCTCACGCCCGGGATGCCCAACCTGCTGTTCCTCTTGGGCGGNNNNNNNNNNGCNGCCNNNGNNNNNNNNNNNNNNGCGCGNGAGNNNN
    +
    !111177777@@@@@@@@@C@@@C@@@CC@@@;@@C@@CC@C@C@@C@C@####!!!!!!!!!!##!###!!!#!!!!!!!!!!!!!!#####!###!!!!
    """

    line_counter = 0
    for line in reads_file:
        line_counter += 1

        # Read every 10th line.
        if line_counter % 40 == 0:
            
            raw_values = map(ord, line.strip())
            normalized_values = map(lambda x: int(x) - 33, raw_values)

            for value in normalized_values:
                #print value
                if value > options.lower_cutoff:
                    quality_values.append(value)
            

            #quality_values.extend(normalized_values)


    #print quality_values
    #print len(quality_values)
    #print '\n'.join(map(str, quality_values))

    summary = tukey_summary(quality_values)
    #print summary

    global cutoff
    cutoff = summary[3]

    # Which function should we use to assign quality?  Good/bad vs. category.
    quality_function = assign_good_bad
    if options.use_category:
        quality_function = assign_category

    decode_function = decode_good_bad
    if options.use_category:
        decode_function = decode_category

    # Go back through sequences and assign good/bad bases.
    reads_file = open(options.reads_filename, 'r')

    new_qualities = []

    for line in reads_file:
        line_counter += 1

        if line_counter % 4 == 0:
            raw_values = map(ord, line.strip())
            normalized_values = map(lambda x: int(x) - options.offset_value, raw_values)
            quals = map(quality_function, normalized_values)
            new_qualities.extend(quals)
            #print quals
            print ''.join(map(decode_function, quals))
        else:
            print line,
            pass

    #print(new_qualities.count("0"))
    #print(new_qualities.count("1"))


def decode_good_bad(value):
    global good_value, offset_value, bad_value
    if value == "1":
        return str(unichr(good_value + offset_value))
    return str(unichr(bad_value + offset_value))


def assign_good_bad(value):
    global cutoff

    if value < cutoff:
        return '0'
    return '1'


def assign_category(value):
    global bin_size

    #print value+1
    #print bin_size
    if value == 0:
        return "1"
    return str(int(math.ceil((value)/bin_size)))

def decode_category(value):
    global bin_size, offset_value

    return str(unichr(int(int(value) * bin_size) + offset_value))

def worst_quality(value):
    # If value is 0, return false.
    if value == 0:
        return False
    return True


def tukey_summary(array):
    """ Given an array of integers, return median, and lower/upper hinge tuple. """

    array.sort()

    median = array[len(array) / 2]
    if not len(array) % 2:
        median = (array[len(array) / 2] + array[len(array) / 2 - 1]) / 2.0

    first_quartile = array[len(array) / 4]
    third_quartile = array[3 * len(array) / 4]

    iqr = third_quartile - first_quartile

    lower_hinge = first_quartile - 1.5 * iqr
    upper_hinge = third_quartile + 1.5 * iqr

    return (median, first_quartile, third_quartile, lower_hinge, upper_hinge)


def get_options():
    parser = OptionParser()
    parser.add_option("-r", "--reads", dest="reads_filename", \
            help="Read file name.")
    parser.add_option("-c", "--cutoff", dest="lower_cutoff", default=2, type="int", \
            help="Don't consider quality values below this value (default = 2).")
    parser.add_option("-g", "--good-value", dest="good_value", default=40, type="int", \
            help="Quality score of a good base.")
    parser.add_option("-b", "--bad-value", dest="bad_value", default=10, type="int", \
            help="Quality score of a bad base.")
    parser.add_option("-o", "--offset-value", dest="offset_value", default=33, type="int", \
            help="Offset score (default 33).")
    parser.add_option("-e", "--encode", dest="encode", default=False, action="store_true", \
            help="Encode.")
    parser.add_option("-d", "--decode", dest="decode", default=False, action="store_true", \
            help="Decode.")
    parser.add_option("-m", "--max_quality_value", dest="max_quality_value", default=40, type="int", \
            help="Max quality value (default = 40).")
    parser.add_option("-n", "--min_quality_value", dest="min_quality_value", default=0, type="int", \
            help="Min quality value (default = 0).")
    parser.add_option("-a", "--categories", dest="num_categories", default=2, type="int", \
            help="Number of categories to use.")
    parser.add_option("-y", "--use-category", dest="use_category", default=False, action="store_true")

    #parser.add_option("-")

    (options, args) = parser.parse_args()

    global offset_value, good_value, bad_value
    offset_value = options.offset_value
    good_value = options.good_value
    bad_value = options.bad_value

    global bin_size
    bin_size = float(options.max_quality_value - options.min_quality_value + 1) / options.num_categories

    return (options,args)


if __name__ == '__main__':
    main()