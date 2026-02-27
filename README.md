# BreakTracer: Tracing inserted sequence fragments at structural variant breakpoints using long-reads

## Installing breaktracer

BreakTracer is available as a [statically linked binary](https://github.com/tobiasrausch/breaktracer/releases/), as a minimal [docker container](https://hub.docker.com/r/trausch/breaktracer/) or as a [singularity containter (SIF file)](https://github.com/tobiasrausch/breaktracer/releases/). The static binaries you can simply download [here](https://github.com/tobiasrausch/breaktracer/releases/) and then make them executable, e.g.:

```
chmod a+x breaktracer_v0.0.5_linux_x86_64bit
./breaktracer_v0.0.5_linux_x86_64bit
```

## Building from source

BreakTracer can be built from source using a recursive clone and make. BreakTracer depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/breaktracer.git`

`cd breaktracer/`

`make all`

## Running breaktracer

BreakTracer has been designed to identify inserted sequence fragments at structural variant (SV) breakpoints using long-read sequencing data. For instance, to identify L1 fragments at SV breakpoints:

`breaktracer find -n L1 -g hg38.fa input.bam > bp.ins.vcf`

BreakTracer can also be used to identify a custom FASTA sequence inserted at SV breakpoints. For instance, to identify a human papillomavirus integration you can use

`breaktracer find -e hpv.seq.fa -g hg38.fa input.bam > bp.ins.vcf`

or with BCF output:

`breaktracer find -e hpv.seq.fa -g hg38.fa -o bp.ins.bcf input.bam`

## License

BreakTracer is free and open source (BSD). Consult the accompanying [LICENSE](https://github.com/tobiasrausch/breaktracer/blob/master/LICENSE) file for more details.


## Credits

[HTSlib](https://github.com/samtools/htslib) is heavily used for alignment processing. [Boost](https://www.boost.org/) for various data structures and algorithms and [Edlib](https://github.com/Martinsos/edlib) for pairwise alignments using edit distance.
