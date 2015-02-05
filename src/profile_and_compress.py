from sklearn.cluster import KMeans
from random import shuffle, sample
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import scipy.misc
import numpy as np
import sys

def qual_ord(phred_string, offset=33):
    return [ord(x)-offset for x in phred_string]

def cluster_rows(matrix):
    linkage_method = 'average' # ['average', 'complete', 'median', 'single', 'ward', 'weighted', 'centroid']:
    D = squareform(pdist(matrix, metric='euclidean'))
    Y = linkage(D, method=linkage_method)
    reorder = dendrogram(Y, no_plot=True)['leaves']
    img = matrix[reorder]
    return img

def sort_rowsum(matrix):
    aa = matrix
    return aa[aa.sum(axis=1).argsort(),:]

N_CLUSTERS = int(sys.argv[3])
TRAINING_SIZE = int(sys.argv[4])

# filename 'gage_rhodo/zzz.quals' 

quals = []
with open(sys.argv[1], 'r') as zzz:
    counter = 1
    for line in zzz:
        if (counter % 4 == 0):
            quals.append(line.rstrip())
        counter += 1

profile_training = sample(quals, TRAINING_SIZE)

features = np.array(map(qual_ord, profile_training))
#print 'starting kmeans'
km = KMeans(n_clusters=N_CLUSTERS, init='random', max_iter=600)
km.fit(features)
#print 'finished kmeans'
if len(sys.argv) > 5:
    scipy.misc.imsave(sys.argv[5], cluster_rows(km.cluster_centers_))

#pred = km.predict(features)
#print ' '.join(map(str,sorted(Counter(pred).values())))
centers = km.cluster_centers_
center_qualstrings = [''.join([chr(max(35,int(x)+33)) for x in center]) for center in centers]


#pred = km.predict(quals)
with open(sys.argv[1], 'r') as zzz:
    with open(sys.argv[2], 'w') as outfile_decomp:
         with open(sys.argv[2]+'.compressed', 'w') as outfile_comp:
             outfile_comp.write('%d\n' % N_CLUSTERS)
             for i_center, center_qualstring in enumerate(center_qualstrings):
                 outfile_comp.write('%d %s\n' % (i_center, center_qualstring))
            
             for i, qq in enumerate(quals):
                 center_class = km.predict(qual_ord(qq))
                 #basis = centers[center_class][0]
                 #qq = np.array(qual_ord(qq))
                 #diff = np.int32(qq - basis)/10
                 #outfile.write(chr(center_class+33) + ''.join(map(chr, diff+20+33)) + '\n')
                 outfile_comp.write(chr(center_class))
                 outfile_decomp.write('%s%s%s%s\n' % (zzz.readline(),zzz.readline(),zzz.readline(),center_qualstrings[center_class]))
                 zzz.readline()

