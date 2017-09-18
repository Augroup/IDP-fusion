IDP-fusion
===

IDP-fusion is an gene Isoform Detection and Prediction tool from Second Generation Sequencing and PacBio sequencing developed by Kin Fai Au. It offers very reliable gene isoform identification with high sensitivity. IDP-fusion can idetify and quantify fusion transcripts.

This is a fork of the original IDP-fusion software.

This manual segment was written for a work in progress set of example data
that didnt quite pan out because I cant fit enough data.  Therefore, I'm going to leave it here as just a hypothetical guide.

You can get the larger test data from the main softwares site.

https://www.healthcare.uiowa.edu/labs/au/IDP-fusion/IDP-fusion_tutorial.asp

If you use this tutorial data you can run the docker as below.  You just need to run the gmap and bowtie aligner building.  Even if one is in there, rebuild it.  It could have been built under a differnet version.  Best to rebuild the indecies.  then modify the run_psl.cfg to point to its missing dependencies.  Then you are ready.

#### 0. Setup the example

First lets get the example.  You should just clone this git repository and use the example data from it. You can clone this repository and we'll work in the example directory.

```bash
$ git clone https://github.com/jason-weirather/IDP-fusion.git
$ cd IDP-fusion/example
$ gunzip data/*.gz
$ ls -lht data
```

#### 1. Correct errors in long reads using short reads

The tools to this are here, but I'll refer you to the IDP for a manual on how to do this. For now lets assume this has been done (and it has) in the example data.

https://github.com/jason-weirather/IDP

#### 2. Align the corrected long reads

! Warning make sure you set `GMAP` to default of multiple paths when you run it for fusion.

Again you can see the `IDP` page if you want to see how do this part.
You could let IDP do this for you, but I caution against it. Its a slow process and the aligners can crash sometimes, so its better to just sort this out now and not deal with it in the IDP run.

We will need a gmap index to use the example

```bash
$ docker run -v $(pwd):/home vacation/idp-fusion \
             gmap_build -d ./gmapindex data/chr20.fa
```

>Byte-coding: 62807695 values < 255, 217826 exceptions >= 255 (0.3%)
>Writing file ./gmapindex/gmapindex.salcpchilddc...done
>Found 3535448 exceptions

Building a gmap index is a slow and memory intesive process on a full genome.  This is just two chromosomes.

We already have these in psl format for deomstration purposes, so you should
skip this step, but you could align them like this (substituting in your corrected lr fasta for the input)

We will use `lr.mm.fa` and `lr.mm.psl` that have already been aligned.

```bash
$ docker run -v $(pwd):/home vacation/idp-fusion \
             gmap --ordered -D ./ -d gmapindex -t 2 -f 1 \
             corrected.lr.fa > corrected.lr.psl
```

#### 3. Align the short reads

I will use `hisat2` to align reads but `runSpliceMap` is included if you want a more classic apporach to the pipeline.

For speed and stability I recommend `hisat2` but it will require an additional processing step on our part.

```bash
$ mkdir hisat2
$ docker run -v $(pwd):/home vacation/idp-fusion \
             hisat2-build data/chr20.fa hisat2/hisat2index
```
>Total time for call to driver() for forward index: 00:01:27

Next we actually align the short reads

```bash
docker run -v $(pwd):/home vacation/idp-fusion \
             hisat2 -x hisat2/hisat2index -U data/sr.fa -f \
| docker run -i vacation/idp-fusion \
             samtools view -Sb - \
| docker run -i -v $(pwd):/home vacation/idp-fusion \
             samtools sort -T midsort - -o sr.sorted.bam
```

>1600000 reads; of these:
>  1600000 (100.00%) were unpaired; of these:
>    48 (0.00%) aligned 0 times
>    1591463 (99.47%) aligned exactly 1 time
>    8489 (0.53%) aligned >1 times
>100.00% overall alignment rate

Looks good!

Unfortunately the IDP component of IDP-fusion needs a different format than the garden variety bam.

To accomodate this we will need to conver the bam into a SpliceMap format sam,
and also create a junction file like SpliceMap does. We use helper scripts for
this part.

First we make the SAM file

```bash
$ docker run -v $(pwd):/home vacation/idp-fusion \
             seq-tools sam_to_splicemap_like_sam sr.sorted.bam \
             -o sr.splicemap-like.sam
```

>processed 1610000 lines. At: chr20:62907177

Next we make the junction file

```bash
$ docker run -v $(pwd):/home vacation/idp-fusion \
             seq-tools sam_to_splicemap_junction_bed sr.sorted.bam \
             -r data/chr20.fa -o sr.splicemap-like.junctions.bed
```

>WARNING skipping non-canonical splice (5/187219)
>finished reading sam

#### 4. Build a bowtie2 index

One of the aligners in here will be looking for a genomic bowtie2 index

```bash
$ mkdir bowtie2/
$ docker run -v $(pwd):/home vacation/idp-fusion \
             bowtie2-build data/chr20.fa bowtie2/bowtie2index
```

#### 5. Run IDP-fusion

On a normal run you will create your own configuration file to describe the run.

Now to actually run IDP-fusion.  This configuration file has been set to use the files created in this example.

In this example we are using an RPKM absolute and fraction cutoff rather than an FDR.  The FDR does not execute well in small datasets or nonmodel organisms.

```bash
$ docker run -v $(pwd):/home vacation/idp-fusion \
             runIDP-fusion.py run.cfg 0
```

Lets see if we have an output .. I didn't wait for all the long reads to correct and align here.

```bash
$ wc -l output/*
```


