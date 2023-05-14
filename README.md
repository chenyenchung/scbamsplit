# scbamsplit

> Last updated: 2023-05-14

An open-source and light-weight tool to subset BAM files based on read tags in parallel
powered by `htslib` and [`uthash`](https://troydhanson.github.io/uthash/) under MIT license.

## Installation

Before installing `scbamsplit`, you will need to install [`samtools`](http://www.htslib.org/) so you get
`htslib` in your environment. You also need [`CMake`](https://cmake.org/) to build the tool.

Once you have both installed, you should run `samtools --version` in the terminal and expect to see:

```
samtools 1.17
Using htslib 1.17
Copyright (C) 2023 Genome Research Ltd.
......
```

Also try `cmake --version`, and you should see this in your terminal:

```
cmake version 3.26.0

CMake suite maintained and supported by Kitware (kitware.com/cmake).
```

If you have both tools ready, get the source code of `scbamsplit` by:

```
git clone https://github.com/chenyenchung/scbamsplit.git
```

or from [this link](https://github.com/chenyenchung/scbamsplit/archive/refs/heads/main.zip).

If you clone the repository, please run the following command in
your terminal:

```bash
cd scbamsplit
mkdir build
cd build
cmake ../ && make
```

If you download from the above link, please decompress it and run the following command in
your terminal:

```bash
cd scbamsplit-man
mkdir build
cd build
cmake ../ && make
```

You should see `scbamsplit` in the build folder, and if you run `./scbamsplit -h`:

```
Program: scbamsplit (Parallel BAM file subset by read tags)
Version: v0.10 (Dependent on htslib v1.17)

Usage: scbamsplit [-f path] [-m path] [-o path] [-q MAPQ] [-t read-tag] [-u read-tag] [-d] [-h] [-n] [-v] 

        -f/--file: the path for input bam file
        -m/--meta: the path for input metadata (an unquoted two-column csv with column names)
        -o/--output: the path to export bam files to (default: ./)
        -q/--mapq: Minimal MAPQ threshold for output (default: 0)
        -t/--filter-tag: The read tag you want to filter against (default: CB)
        -u/--umi-tag: The read tag containing UMIs (default: UB)
        -d/--dedup: Remove duplicated reads with the same cell barcode/UMI combination
        -n/--dry-run: Only print out parameters
        -v/--verbose: Print out parameters
        -h/--help: Show this documentation


```

This executable can be moved or soft linked to a directory in your path (e.g., `$HOME/.local/bin`) so you may
invoke it anywhere as a command.

## Usage

`scbamsplit` requires:

1. a BAM file with reads that are marked by read tags (-f/--file)
2. a comma-separated file (.csv) in which the first column is the value of the tag to filtered
while the second is the subset identity (e.g., cluster, sample...). `scbamsplit` will generate a BAM file
for each identity and export reads that contains a tag that belongs to this identity in the file.
3. An output directory (-o/--output) to put the exported BAM files.

By default, the read tag used to query the provided metadata is `CB` (the read tag that contains
corrected cell barcode in `cellranger`-aligned BAMs). To use other read tags, set `-t` (e.g., if
tag `SM` contains sample information that you want to use to split the file, try `-t SM` or
`--filter-tag SM`).

It is advisable to add `-n`/`--dry-run` to the command, so when you execute the command, nothing but
**printing the information for the run** will happen. If you are happy with what you see, remove
`-n`/`--dry-run` to really run the tool.

`scbamsplit` is pretty easy on memory usage and should only use less than 1GB (size of metadata +
iterating through the input BAM with `htslib`) in most use cases. However, if you enabled
the optional deduplication feature (please see below), memory usage will increase linearly with
your library complexity. For a regular 10X Genomics gene expression run with 10k cells, I would
recommend running with 8 - 16 GB to begin with and adjust accordingly based on the test run.


There are some extra functionalities that are optional:

### MAPQ filtering

If you want to set an lower limit of MAPQ for reads to be exported, try `-q`
(e.g., if you want only reads with MAPQ>=30, try `-q 30` or `--mapq 30`).

### UMI-based deduplication

Some sequencing techniques involves adding a unique molecule index (UMI) to
each fragments of interest. In such cases, one might want to keep one read
from each unique fragment once. To enable this function, please use
`-d` (or `--dedup`) flag to turn on deduplication.

When deduplication is on, the combinations of **filter tag value (e.g., cell barcode)** and
**UMIs** will be stored, and only the first read with the same combination will be exported.

The default read tag for UMIs is `UB`, and you may change this with `-u` or `--umi-tag` (
e.g., to set UMI tag to `UM`, try `-u UM` or `--umi-tag UM`).

Please note that this function assumes UMIs and filter tags are corrected and contains
no PCR errors.
