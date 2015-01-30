from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from scipy.misc import toimage, imread
import numpy as np
import StringIO
from PIL import Image

def qual_ord(phred_string, offset=33):
    return [ord(x)-offset for x in phred_string]

def png_to_reads(png_file, order_file):
    orders = []
    with open(order_file) as in_order:
        for line in in_order:
            orders.append(map(int,line.rstrip().split()))
    image = Image.open(png_file)
    quals = np.asarray(image)
    original_ordering = np.zeros(quals.shape, dtype=np.int8)
    curr_pos = 0
    for current_block in orders:
        for i, o in enumerate(current_block):
            original_ordering[curr_pos+o, :] = quals[curr_pos+i, :]
        curr_pos += len(current_block)
    return original_ordering, image.getpalette()[::3]

def matrix_to_qualstrings(matrix, palette, offset=33):
    def convert_with_palette(number):
        return chr(palette[number]+offset)

    matrix = matrix
    for row in matrix:
        print ''.join(map(convert_with_palette, row))

for i in range(400):
    #mx, pal = png_to_reads(('clustered_pngs/wed_256/compressed/wed_256_%d-or8.png' % i), ('clustered_pngs/wed_256/wed_256_%d.txt' % i))
    mx, pal = png_to_reads(('clustered_pngs/sj_2/sj_2_%d-or8.png' % i), ('clustered_pngs/sj_2/sj_2_%d.txt' % i))
    matrix_to_qualstrings(mx, pal)