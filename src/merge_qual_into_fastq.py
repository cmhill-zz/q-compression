import sys

with open(sys.argv[1], 'r') as in_fastq:
    with open(sys.argv[2], 'r') as in_qual:
        while True:
            try:
                line1 = next(in_fastq)
                line2 = next(in_fastq)
                line3 = next(in_fastq)
                _ = next(in_fastq)
                line4 = next(in_qual)
                sys.stdout.write(line1 + line2 + line3 + line4)
            except StopIteration:
                break
