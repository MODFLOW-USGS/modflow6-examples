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
- [Contributing examples](#contributing-examples)
- [Releasing the examples](#releasing-the-examples)
  - [Extract notebooks and tables from scripts](#extract-notebooks-and-tables-from-scripts)
  - [Build the documentation](#build-the-documentation)

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

The examples can be run with Jupyter or with Pytest.

### Using `jupyter`

To start a Jupyter browser interface, run `jupyter notebook` from the `notebooks/` directory.

### Using `pytest`

Pytest can be used to run all of the examples from the `scripts` directory. For instance, to run the examples on as many cores as your machine can spare, issue the following command:

```commandline
pytest -v -n auto
```

Remove `-n auto` to run in serial instead of parallel.

## Contributing examples

Adding a new example requires adding a new example script in the `scripts/` folder and adding a new LaTeX problem description in the `doc/sections/` folder. After a new example is added, a new release should be made (see below).

## Releasing the examples

A new release is automatically created whenever code is merged into the trunk branch of this repository. Steps to manually prepare for a release are listed below.

### Extract notebooks and tables from scripts

The example scripts are converted into jupyter notebooks, and latex tables are created from the scripts using jupytext.  To convert all of the scripts, run the following command in the scripts directory:

```commandline
python process-scripts.py
```

### Build the documentation

If the figures and tables were generated correctly, then it should be possible to build the pdf document for the examples. The pdf is created by converting the `doc/mf6examples.tex` document into a pdf with `pdftex` or a similar tool.
