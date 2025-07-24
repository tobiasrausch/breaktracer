# BreakTracer: Tracing inserted sequence fragments at structural variant breakpoints

## Installing breaktracer

BreakTracer can be build from source using a recursive clone and make. BreakTracer depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/breaktracer.git`

`cd breaktracer/`

`make all`

## Running breaktracer

To identify L1 fragments at structural variant breakpoints.

`breaktracer find -g hg38.fa input.bam > bp.reads.tsv`

License
-------
BreakTracer is free and open source (BSD). Consult the accompanying [LICENSE](https://github.com/tobiasrausch/breaktracer/blob/master/LICENSE) file for more details.


Credits
-------
[HTSlib](https://github.com/samtools/htslib) is heavily used for alignment processing. [Boost](https://www.boost.org/) for various data structures and algorithms and [Edlib](https://github.com/Martinsos/edlib) for pairwise alignments using edit distance.
