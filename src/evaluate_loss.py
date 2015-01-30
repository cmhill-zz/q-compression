import numpy as np
from collections import Counter
import sys
from itertools import izip

# awk 'NR%4==0' whatever.fastq > quals_only.qual

def qual_ord(phred_string, offset=33):
    return map(lambda x: (ord(x)-offset), phred_string)


def evaluate_loss(qstring1, qstring2):
    q1, q2 = np.array(qual_ord(qstring1)), np.array(qual_ord(qstring2))
    diff = q1 - q2
    mse = (diff ** 2).mean()
    return mse, Counter(q1 - q2)

file1 = sys.argv[1]
file2 = sys.argv[2]

mses = []
diff_histogram = Counter([])

IS_FASTQ = False

j = 0

with open(file1, 'r') as in1:
    with open(file2, 'r') as in2:
        for i, (line1, line2) in enumerate(izip(in1, in2)):
            if IS_FASTQ and (i%4 != 3):
                continue
            # if (j%10000) == 0:
            #     print j
            j += 1
            if (j % 50) != 0:  # only sample every 50th read. Otherwise quite slow
                continue
            mse, diffcounts = evaluate_loss(line1.rstrip(), line2.rstrip())
            mses.append(mse)
            diff_histogram += diffcounts

mses = np.array(mses)
print 'mean/std mse', mses.mean(), mses.std()
#hist = np.histogram(mses, bins=50)
#print 'error histogram', diff_histogram

