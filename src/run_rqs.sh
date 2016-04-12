#!/bin/bash
# echo First convert FASTQ files to a fake SAM format:
$RQS/fastq2fsam.sh $1
output=$2/`basename $1`.fsam
mv $1.fsam $output

# echo Generate our dictionary with read multiplicity r=2:
#../generate_dict 2 dict.db $1.fsam
$RQS/generate_dict 2 $2/dict.db $output

# echo In this example, we compress the corpus we used to generate the dictionary:
$RQS/sparsify $2/dict.db $output

# echo After sparsification, we still need to cut-off high quality scores:
# echo "Let's choose a threshold of 'I', or Phred quality 40 under most encodings"
# echo "../threshold 'I' *.fsam.filtered"
$RQS/threshold 'I' $output.filtered

# echo Compression is now done, but let\'s see how well we did using BZIP2:
# echo "We'll cut out column 11 (the quality scores) and pipe it through BZIP2"
# echo and then compute the bits per quality value needed.
for filename in $output.filtered.reduced
do
        file=`basename $filename`
        cut -f11 $filename > $2/$file.qual
        orig_size=`wc -c < $2/$file.qual`
        orig_lines=`wc -l < $2/$file.qual`
        orig_size=`echo "$orig_size - $orig_lines" | bc`
        bzip2 -f $2/$file.qual
        new_size=`wc -c < $2/$file.qual.bz2`
        #rm $file.qual.bz2
        echo -e $file:'\t' `echo "scale=4; 1/( $orig_size / ( $new_size * 8)) " | bc` bits / quality score
done

# echo "Although we're basically done, you might want to convert your files back"
# echo "from this fake SAM format to a FASTQ file:"
# echo "../fsam2fastq.sh *.filtered.reduced"
$RQS/fsam2fastq.sh $output.filtered.reduced
mv $output.filtered.reduced.fastq $2/`basename $1`
