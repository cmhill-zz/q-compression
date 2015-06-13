Lossy compression of quality values
=============
A pipeline for evaluating several different lossy compression methods for quality values, and the resulting impacts on common bioinformatic analyses using sequence read data: preprocessing, genome assembly, and alignment of short reads to a reference sequence.

Dependencies
-------
* Python 2.7.6
* R 3.1.2
* Optional:
  * [Sickle](https://github.com/najoshi/sickle)
  * [FastQ-tools](https://github.com/dcjones/fastq-tools)
  * [QualComp](http://web.stanford.edu/~iochoa/QualComp.html)
  * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [Allpaths-LG](https://www.broadinstitute.org/software/allpaths-lg/blog/?page_id=12)
  * [LAP](http://assembly-eval.sourceforge.net/manual.shtml)

Usage
------
```
  Usage: compress.py [options]

  Options:
    -h, --help            show this help message and exit
    -r UNPAIRED_READS_FILENAMES, --reads=UNPAIRED_READS_FILENAMES
                          Unpaired FASTQ filenames.
    -1 FIRST_MATE_FILENAMES, --first=FIRST_MATE_FILENAMES
                          First mate FASTQ filenames.
    -2 SECOND_MATE_FILENAMES, --second=SECOND_MATE_FILENAMES
                          Second mate FASTQ filenames.
    -f REFERENCE_FASTA, --reference=REFERENCE_FASTA
                          Reference FASTA filename.
    -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                          Output directory.
    -a, --assemble        Run assembly evaluation
    -p, --preprocessing   Run preprocessing tools evaluation
    -b, --alignment       Run alignment evaluation (using Bowtie2).
    -s, --sort            Sort FASTQ reads before the pipeline begins (requires
                          fastq-sort).
    --poly-degrees=POLY_DEGREES
                          Comma-separated list of polynomial degrees to use for
                          regression.
    --training-size=TRAINING_SIZE
                          Training size used for clustering.
    --profile-sizes=PROFILE_SIZES
                          Comma-separated list of number of profiles to use.
    --rates=RATES         QualComp parameter for setting the  bits/reads.
    --clusters=CLUSTERS   QualComp parameter for setting number of clusters.
    --max-qv=MAX_QUALITY  Use this value for max quality value compression.
    --min-qv=MIN_QUALITY  Use this value for min quality value compression.
    -t THREADS, --threads=THREADS
                          Number of threads (default 32
```

Getting started
------

Clone the repository and run the pipeline on the included test data set:
```
python ./src/compress.py -r test/frag_small_1.fastq  -o test_results/ --training-size 100 --profile-sizes 64 --poly-degrees 0,1 -t 4 -s --rate 6,10,100
```

This will compress the frag_small_1.fastq file in the following ways:
* Run profile-based compression using [--profiles-sizes] profiles calculated from a training set of size [--training-size].
* Run polynomial regression compression using [--poly-degrees]-order polynomials.
* Run QualComp with rates [--rate].
* Goodbad encoding is run by default.

The results will be found in the directory
* test_results/results/frag_small_1.fastq.compression

```
compression     bases   orig_size       orig_bzip2      comp_size       comp_bzip2      mse     bits_bp
original        303000  306000  82728   NA      NA      0.0     2.18423762376
goodbad 303000  306000  12102   37876   14471   33.5415841584   0.382072607261
maxqual 303000  306000  51      37876   14471   666.901320132   0.382072607261
minqual 303000  306000  51      37876   14471   702.148844884   0.382072607261
degree_0        303000  306000  2222    24000   2195    151.408580858   0.0579537953795
degree_1        303000  306000  14206   48000   32786   47.8793729373   0.865636963696
profile_64      303000  306000  7587    9464    5614    21.4202970297   0.148224422442
qualcomp_r6     303000  306000  35583   2272    2678    20.252970297    0.0707062706271
```

Additional evaluations
------
In addition to evaluating the compression of each strategy, we provide pipelines to evaluate the effect of lossy compression on **read alignment** (using Bowtie2), **preprocessing** (using Sickle) and **assembly** (using Allpaths-LG).

### Read alignment
Include the ```-b``` and ```-f [REFERENCE.FASTA]``` parameter. ```bowtie2``` and ```bowtie2-build``` must be found in the ```$PATH``` environmental variable. Results will be located in ```[OUTPUT_DIR]/results/[READS].alignment```.  The compressed sequences are aligned to the reference genome and compared with the alignments of the original sequences.  For each compression method we store the total amount of mapped sequences (```mapped```), the number of sequences aligned by the compression strategy and original (```shared```), the alignments unique to the original (```uniq_orig```), the alignments unique to the compression strategy (```uniq_comp```), and lastly, the number of sequences unaligned by both methods (```unmapped_shared```).

```
compression     mapped  shared  uniq_orig       uniq_comp       unmapped_shared
original        892613  892613  0       0       132821
goodbad 892050  891864  749     186     132635
maxqual 746716  746716  145897  0       132821
minqual 903434  892613  0       10821   122000
degree_0        851749  851682  40931   67      132754
degree_1        883445  883390  9223    55      132766
degree_3        889654  889537  3076    117     132704
degree_5        891174  891019  1594    155     132666
degree_7        891633  891479  1134    154     132667
profile_64      891897  891753  860     144     132677
profile_128     892095  891952  661     143     132678
profile_256     892170  892051  562     119     132702
qualcomp_r6     891679  891375  1238    304     132517

```

### Preprocessing
Include the ```-p``` option. ```sickle``` must be found in ```$PATH```. Results will be located in ```[OUTPUT_DIR]/results/[READS].preprocessing```. The results contain information comparing the results of running the preprocessing tool ```sickle``` on the original sequences with the compressed version. For each compression strategy, we store the sequences only kept in the original (```uniq_orig```), sequences only kept by the compressed strategy (```uniq_comp```), sequences kept by both methods (```common_kept```), and sequences discarded by both methods (```common_discard```).  Lastly, the total amount of bases kept by the compression strategy is displayed (```bases_kept```).


```
compression     uniq_orig       uniq_comp       common_kept     common_discard  bases_kept
original        0       0       2561    439     187894
goodbad 53      6       2508    433     183606
maxqual 0       439     2561    0       303000
minqual 2561    1       0       438     0
degree_0        514     2       2047    437     206949
degree_1        38      35      2523    404     188650
profile_64      1       63      2560    376     192918
qualcomp_r6     13      45      2548    394     190482
qualcomp_r10    13      44      2548    395     191504
```

### Assembly
If the supplied FASTQ files fit the Allpaths-LG naming scheme (frag_[12].fastq, shortjump_[12].fastq, etc.), then include ```-a``` parameter to perform lossy compression and assemble the compressed reads using Allpaths-LG.  ```-f [REFERENCE.FASTA]``` must also be included to perform reference-based comparisons. ***Important:*** the analysis scripts must be downloaded from the [GAGE website](http://gage.cbcb.umd.edu/results/) and the directory must be saved in an environmental variable named ```MUMMER```.

```
compression     ref_bases       asm_bases       N50     missing_ref_bases       avg_idy snps    indels_gt5bp    inversions      reloc   transloc        corr_N50        LAP
original        4603060 4607865 3192325 25229   99.98   222     170     0       9       2       40171   -34.9185434849
goodbad 4603060 4625295 3171102 35245   99.96   179     283     0       11      2       22973   -34.9795497955
maxqual 4603060 4761628 182153  682756  99.95   730     612     2       12      7       4976    -39.7830873617
minqual 4603060 607423  0       4097747 99.94   62      78      2       10      15      0       -68.8205252116
degree_0        4603060 4637494 911256  209868  99.95   655     451     2       14      7       11727   -36.2046363265
degree_1        4603060 4617120 3195723 66808   99.95   541     363     2       12      4       17052   -35.4272904581
degree_3        4603060 4612062 3170886 29715   99.97   302     217     0       6       3       31376   -35.0250248098
```
