# LongRNAqc

## Read artifacts and isoform matching

(LongRNAqc)[] is a workflow intended to look at base reads quality rather than collapsed reads, and to easily compare multiple samples on Terra. Taking an aligned BAM as input for a sample, it optionally runs SQANTI3 (which normally does not have the option to use existing alignments) with some parallelization, LRAA and IsoQuant.

Once multiple samples have been processed, a (downstream plotting workflow)[] can be run on a selection of samples to generate a comparison report based on the outputs. Currently this only works based on Sqanti3 outputs, but the outputs of LRAA are extended, it will be supported as well.


## Sequencing quality

(A standalone workflow)[] taking an aligned BAM is also available to look at error rates and profiles between samples. However to be able to produce stats about mismatch errors, the CIGAR string in the alignments need to use the "=" (match) and "X" (mismatch) operators rather than the "M" (match or mistmatch) operator. With Minimap2, this is done by including the "--eqx" argument.



