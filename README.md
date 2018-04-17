[![Build Status](https://travis-ci.org/Cibiv/gwpcR.svg?branch=master)](https://travis-ci.org/Cibiv/gwpcR)
[![Coverage Status](https://coveralls.io/repos/github/Cibiv/gwpcR/badge.svg?branch=master)](https://coveralls.io/github/Cibiv/gwpcR?branch=master)
[![Install with Conda](https://anaconda.org/bioconda/r-gwpcr/badges/installer/conda.svg)](https://anaconda.org/bioconda/r-gwpcr)

# Description

__Motivation__: Counting molecules using next-generation sequencing (NGS) suffers from PCR amplification bias, which reduces the accuracy of many quantitative NGS-based experimental methods such as RNA-Seq. This is true even if molecules are made distinguishable using unique molecular identifiers (UMIs) before PCR amplification, and distinct UMIs are counted instead of reads: Molecules that are lost entirely during the sequencing process will still cause under-estimation of the molecule count, and amplification artifacts like PCR chimeras create phantom UMIs and thus cause over-estimation.

__Results__: *gwpcR* implements mechanistic model of PCR amplification that allows correction of both types of errors. In our [paper](https://www.biorxiv.org/content/early/2017/11/13/217778) we demonstrate that the model describes UMI-based NGS experiments well, and that using it to filter phantoms and correct for lost molecules considerably increases the accuracy of measured molecule counts over just counting the number of distinct UMIs.

__Using *gwpcR*__: The easiest way to integrate our loss- and phantom-correction algorithm into your UMI pipeline is by using our command-line tool [TRUmiCount](https://cibiv.github.io/trumicount) (based on `gwpcR` of course). [TRUmiCount](https://cibiv.github.io/trumicount) integrates with [UMI-Tools](https://github.com/CGATOxford/UMI-tools), and allows you to get from a BAM file containing mapped reads to a per-gene count table already corrected for sequencing errors, amplification artifacts and lost molecules with a single command.

__Using *gwpcR*__ directly: If your pipeline is already R-based, you might want to integrate *gwpcR* directly instead of using our command-line tool [TRUmiCount](https://cibiv.github.io/trumicount). After installing the *gwpcR* package, see `help(gwpcrpois.est)` and `help(gwpcrpois.groupest)`.

# Installation

## Using Conda 

If you're already using [Conda](https://bioconda.github.io/), you can install TRUmiCount from the [Bioconda](https://conda.io/) channel by doing
```
  conda install -c bioconda r-gwpcr
```
## Using *devtools*

*gwpcR* can also be installed directly from this GitHub repository using R's *devtools* package.
```
  install.packages("devtools")
  devtools::install_github("Cibiv/gwpcR", ref="latest-release")
```

# Publications

The implemented model is described in detail in our paper

Florian G. Pflug, Arndt von Haeseler. (2018). TRUmiCount: Correctly counting
absolute numbers of molecules using unique molecular identifiers.
<i>Bioinformatics</i>, DOI: https://doi.org/10.1093/bioinformatics/bty283.

# License

*gwpcR* is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

*gwpcR* is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more
details.
