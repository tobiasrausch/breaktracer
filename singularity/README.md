You can build a [breaktracer](https://github.com/tobiasrausch/breaktracer) singularity container (SIF file) using

`sudo singularity build breaktracer.sif breaktracer.def`

Once you have built the container you can run analysis using

`singularity exec breaktracer.sif breaktracer --help`
