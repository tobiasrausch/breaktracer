# BreakTracer: Tracing inserted sequence fragments at structural variant breakpoints using long-reads

## Installing breaktracer

BreakTracer is available as a [statically linked binary](https://github.com/tobiasrausch/breaktracer/releases/), as a minimal [docker container](https://hub.docker.com/r/trausch/breaktracer/) or as a [singularity containter (SIF file)](https://github.com/tobiasrausch/breaktracer/releases/).

```
wget https://github.com/tobiasrausch/breaktracer/releases/download/v0.0.3/breaktracer_v0.0.3_linux_x86_64bit
chmod a+x breaktracer_v0.0.3_linux_x86_64bit
./breaktracer_v0.0.3_linux_x86_64bit
```

## Building from source

BreakTracer can be built from source using a recursive clone and make. BreakTracer depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/breaktracer.git`

`cd breaktracer/`

`make all`

## Running breaktracer

BreakTracer has been designed to identify inserted sequence fragments at structural variant (SV) breakpoints using long-read sequencing data. For instance, to identify L1 fragments at SV breakpoints:

`breaktracer find -t L1 -g hg38.fa input.bam > breakpoint.insertion.tsv`

BreakTracer can also be used to identify a custom FASTA sequence inserted at SV breakpoints. For instance, to identify a human papillomavirus integration:

`breaktracer find -e hpv.seq.fa -g hg38.fa input.bam > breakpoint.insertion.tsv`

## Building a genome mask

Sequences that are inserted but also present in the reference genome can generate complex mapping artifacts (depending on the aligner). To account for these mapping artifacts, please create a genome mask for such inserted sequences, e.g., for the standard L1 sequence:

`breaktracer mask hg38.fa > hg38.L1.mask.fa`

`samtools faidx hg38.L1.mask.fa`

This genome mask only needs to be created once for a given reference genome. The search uses the alignment edit distance and therefore takes approximately ~10 hours for a human genome.

## Running breaktracer with a genome mask

To remove SV breakpoints that correspond to a region in the reference genome that matches the inserted sequence, use a pre-built genome mask.

`breaktracer find -t L1 -g hg38.fa -k hg38.L1.mask.fa input.bam > breakpoint.insertion.tsv`

## License

BreakTracer is free and open source (BSD). Consult the accompanying [LICENSE](https://github.com/tobiasrausch/breaktracer/blob/master/LICENSE) file for more details.


## Credits

[HTSlib](https://github.com/samtools/htslib) is heavily used for alignment processing. [Boost](https://www.boost.org/) for various data structures and algorithms and [Edlib](https://github.com/Martinsos/edlib) for pairwise alignments using edit distance.
