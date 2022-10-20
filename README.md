This is the repository for `Sgootr`, a tool that jointly infers a tumor lineage tree and selects lineage-informative features using single-cell methylation sequencing data.

![Schema Figure for Sgootr](/assets/sysarch.png)


# Getting Started

Follow [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install `Snakemake`.
Follow [instructions](https://www.r-project.org/) to install `R`.

```console
-$ git clone https://github.com/liuy0421/Sgootr.git
-$ cd Sgootr
```

Download `data.tar.gz` into `Sgootr/` from [here](https://umd.box.com/v/sgootr-crc01).

```console
-$ tar -xf data.tar.gz
-$ snakemake --cores <number of cores> --use-conda
```
