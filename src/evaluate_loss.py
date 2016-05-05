import math
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

    # L1
    l1 = np.linalg.norm((q1 - q2), ord=1)

    # Lorentzian
    lorentzian = math.log(1 + l1, 2)

    return mse, l1, lorentzian


file1 = sys.argv[1]
file2 = sys.argv[2]

mses = []
L1s = []
lorentzians = []

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

            mse, L1, lorentzian = evaluate_loss(line1.rstrip(), line2.rstrip())
            mses.append(mse)
            L1s.append(L1)
            lorentzians.append(lorentzian)
            #diff_histogram += diffcounts

mses = np.array(mses)
L1s = np.array(L1s)
lorentzians = np.array(lorentzians)

print mses.mean(),'\t',mses.std(),'\t',L1s.mean(),'\t',L1s.std(),'\t',lorentzians.mean(),'\t',lorentzians.std()
#hist = np.histogram(mses, bins=50)
#print 'error histogram', diff_histogram
