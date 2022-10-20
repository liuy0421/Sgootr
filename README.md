This is the repository for `Sgootr`, a tool that jointly infers a tumor lineage tree and selects lineage-informative features using single-cell methylation sequencing data.

![Schema Figure for Sgootr](/assets/sysarch.png)


# Getting Started

Follow [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install `Snakemake`.
Follow [instructions](https://www.r-project.org/) to install `R`.
Then:

```console
-$ git clone https://github.com/liuy0421/Sgootr.git
-$ cd Sgootr
```

To apply `Sgootr` to multiregionally-sampled metastatic colorectal patient CRC01 made available by Bian *et al.* [^1], download `data.tar.gz` into `Sgootr/` from [here](https://umd.box.com/v/sgootr-crc01), then:
[^1] Bian, S., Hou, Y., Zhou, X., Li, X., Yong, J., Wang, Y., Wang, W., Yan, J., Hu, B., Guo, H., Wang, J.,
Gao, S., Mao, Y., Dong, J., Zhu, P., Xiu, D., Yan, L., Wen, L., Qiao, J., Tang, F., Fu, W.: Single-cell multiomics sequencing and analyses of human colorectal cancer. Science **362**(6418), 1060-1063 (Nov 2018). [https://doi.org/10.1126/science.aao3791](https://doi.org/10.1126/science.aao3791)

```console
-$ tar -xf data.tar.gz
-$ snakemake --cores <number of cores> --use-conda
```
