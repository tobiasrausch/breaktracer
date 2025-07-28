# BreakTracer: Tracing inserted sequence fragments at structural variant breakpoints

## Installing breaktracer

BreakTracer is available as a [statically linked binary](https://github.com/tobiasrausch/breaktracer/releases/) or as a [singularity containter (SIF file)](https://github.com/tobiasrausch/breaktracer/releases/).

`wget https://github.com/tobiasrausch/breaktracer/releases/download/v0.0.3/breaktracer_v0.0.3_linux_x86_64bit`

`chmod a+x breaktracer_v0.0.3_linux_x86_64bit`

`./breaktracer_v0.0.3_linux_x86_64bit`

## Building from source

BreakTracer can be built from source using a recursive clone and make. BreakTracer depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/breaktracer.git`

`cd breaktracer/`

`make all`

## Running breaktracer

To identify L1 insertion fragments at structural variant breakpoints.

`breaktracer find -t L1 -g hg38.fa input.bam > breakpoint.insertion.tsv`

To identify a custom FASTA sequence inserted at structural variant breakpoints.

`breaktracer find -e insertion.seq.fa -g hg38.fa input.bam > breakpoint.insertion.tsv`

## License

BreakTracer is free and open source (BSD). Consult the accompanying [LICENSE](https://github.com/tobiasrausch/breaktracer/blob/master/LICENSE) file for more details.


## Credits

[HTSlib](https://github.com/samtools/htslib) is heavily used for alignment processing. [Boost](https://www.boost.org/) for various data structures and algorithms and [Edlib](https://github.com/Martinsos/edlib) for pairwise alignments using edit distance.
