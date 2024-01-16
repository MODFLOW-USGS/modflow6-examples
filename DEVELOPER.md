# Developing MODFLOW 6 examples

This document describes how to set up an environment to use or develop the MODFLOW 6 example models.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Overview](#overview)
- [Setup](#setup)
  - [Install MODFLOW 6 and related programs](#install-modflow-6-and-related-programs)
  - [Install Python dependencies](#install-python-dependencies)
  - [Update FloPy](#update-flopy)
- [Running the examples](#running-the-examples)
  - [Using `jupyter`](#using-jupyter)
  - [Using `pytest`](#using-pytest)
  - [Using `jupytext`](#using-jupytext)
- [Contributing examples](#contributing-examples)
- [Releasing the examples](#releasing-the-examples)
  - [Generate notebooks and tables](#generate-notebooks-and-tables)
  - [Build examples documentation](#build-examples-documentation)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Overview

The examples are Python scripts managed with [`jupytext`](https://github.com/mwouts/jupytext). The Python files in `scripts/` are the ultimate source of truth. Data and notebooks are generated from the scripts. There are more examples than scripts: some of the scripts produce multiple examples.

## Setup

This section lists steps to set up a development environment, build the examples from scripts, and finally build the examples documentation.

- install MODFLOW 6 and related programs
- Install Python dependencies
- make sure FloPy is up to date

### Install MODFLOW 6 and related programs

Install modeling programs using the `get-modflow` utility that is available from FloPy. From any directory, issue the following commands:

```commandline
get-modflow :
get-modflow --repo modflow6-nightly-build :
```

The command will show a prompt to select the install location. Any location can be selected, as `get-modflow` will put the executables on your path (the notebooks expect binaries to be available on the path).

### Install Python dependencies

Several Python packages are required to create the MODFLOW 6 examples. These are listed in `etc/requirements.pip.txt` and `etc/requirements.usgs.txt`.

### Update FloPy

It's important that FloPy is up-to-date with the latest changes to MODFLOW 6. The latest release of FloPy is always up to date with the latest release of MODFLOW 6.

To manually update FloPy from some branch of the [official MODFLOW 6 repository](https://github.com/MODFLOW-USGS/modflow6):

```commandline
import flopy
flopy.mf6.utils.generate_classes(owner="MODFLOW-USGS", branch="master", backup=True)
```

The above is equivalent to calling the function with no arguments. Arguments may be substituted to select an alternative repository (e.g. your fork) or branch.

## Running the examples

The example scripts generate input and output files in subdirectories of `examples/`. Three environment variables control their behavior, each only effective if the previous are, all accepting case-insensitive `true/false`:

- `RUN`: whether to run models &mdash; set `false` to build but not run
- `PLOT`: whether to plot model results
- `SAVE`: whether to save plots to file

The example scripts can be run directly, tested via `pytest`, converted to notebooks and/or executed with `jupytext`, or run as notebooks with `jupyter`. When the scripts are run directly, all three configuration variables default to true. When run via `pytest`, all three default to false. There are corresponding `pytest` CLI arguments `--run`, `--plot`, and `--save`.

For instance, to test the example scripts in parallel, running models but omitting plots:

```shell
pytest -n auto test_scripts.py --run
```

To test notebook conversion/execution with `jupytext`:

```shell
pytest -n auto test_notebooks.py --run
```

To convert an example script to a notebook manually with `jupytext`:

```shell
jupytext --from py --to ipynb scripts/ex-gwf-twri.py -o notebooks/ex-gwf-twri.ipynb --execute
```

To start a Jupyter browser interface, run `jupyter notebook` from the `notebooks/` directory.

## Contributing examples

Adding a new example requires:

* adding a new example script in the `scripts/` folder
* adding a new LaTeX problem description in the `doc/sections/` folder.

Then open a pull request from your fork of the repository.

## Releasing the examples

A new release is automatically created whenever code is merged into the trunk branch of this repository. Steps to manually prepare for a release are listed below.

1. Generate notebooks and tables
2. Build examples documentation

### Generate notebooks and tables

The example scripts must be converted into jupyter notebooks using `jupytext`, then latex tables generated from specially formatted code/comments. Run `process-scripts.py` in the `scripts/` directory to do this.

### Build examples documentation

If the figures and tables were generated correctly, then it should be possible to build the examples PDF document. The PDF can be created by processing `doc/mf6examples.tex` with `pdftex` or a similar tool.
