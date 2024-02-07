# Developing MODFLOW 6 examples

This document describes development procedures and conventions for the MODFLOW 6 example models.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Overview](#overview)
- [Prerequisites](#prerequisites)
  - [Create a Python environment](#create-a-python-environment)
  - [Install MODFLOW programs](#install-modflow-programs)
  - [Update FloPy classes](#update-flopy-classes)
- [Running the examples](#running-the-examples)
  - [Running with `pytest`](#running-with-pytest)
  - [Running with `jupytext`](#running-with-jupytext)
  - [Running with `jupyter`](#running-with-jupyter)
- [Contributing examples](#contributing-examples)
- [Releasing the examples](#releasing-the-examples)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Overview

The example models in this repository are composed with FloPy in the Python files under `scripts/`. Notebooks are generated from these with [`jupytext`](https://github.com/mwouts/jupytext) for the [online documentation](https://modflow6-examples.readthedocs.io/en/stable/) and/or for use with `jupyter`. A small volume of supporting data live in files under `data/`.

**Note:** example scripts do not access the `data/` folder directly &mdash; rather, on first run, they download and cache files from GitHub with a tool called [Pooch](https://www.fatiando.org/pooch/latest/). This allows running scripts and/or notebooks downloaded from the [documentation site](https://modflow6-examples.readthedocs.io/en/stable/) without first cloning the repository and/or manually downloading data files.

Scripts may contain one or more example scenarios. Each script creates one or more subdirectories of `examples/` to use as a workspace, one per example scenario. Internally, the scripts are structured similarly, proceeding through the following steps:

1. Define & build models
2. Write input files
3. Run models
4. Plot results

The scripts parse environment variables to control several behaviors, all of which accept case-insensitive `true/false` strings:

- `WRITE`: whether to write model input files to subdirectories of `examples/`
- `RUN`: whether to run models
- `PLOT`: whether to create plots of results
- `PLOT_SHOW`: whether to show static plots (disabled by default when run with `pytest`)
- `PLOT_SAVE`: whether to save static plots to `figures/`
- `GIF`: whether to save GIFs to `figures/` (only relevant for a small subset of scripts)

If variables are not found when a script is run directly, behaviors are enabled by default. When scripts are run via `pytest`, by default plots are not shown (to avoid the need to manually close plot widgets). Pytest CLI flags can be used to control each switch &mdash; see below for details.

## Prerequisites

To develop the examples one must first:

- create a Python environment
- install MODFLOW programs
- update FloPy classes

### Create a Python environment

Several Python packages are required to run the example scripts. These are listed in `etc/requirements.pip.txt` and `etc/requirements.usgs.txt`. Once a Python environment has been created, e.g. with `venv` or Conda, dependencies can be installed with:

```shell
pip install -r etc/requirements.pip.txt
pip install -r etc/requirements.usgs.txt
```

### Install MODFLOW programs

Besides MODFLOW 6, the examples use a number of MODFLOW-related programs, discovered on the system `PATH`.

With FloPy installed in the active Python environment, these can all be installed with:

```shell
get-modflow :
```

You will be prompted to select an install location &mdash; your virtual environment's bindir is a good place, as it is already on the system `PATH` and keeps things isolated.

Typically one develops against the cutting edge of MF6 &mdash; this might entail adding a local development version of MF6 to the bindir as well. Alternatively, the MF6 nightly build can be installed with:

```shell
get-modflow --repo modflow6-nightly-build :
```

See the [FloPy documentation](https://flopy.readthedocs.io/en/stable/md/get_modflow.html) for more info on the `get-modflow` utility.

### Update FloPy classes

FloPy and MODFLOW 6 versions must be kept in sync for FloPy to properly generate and consume MF6 input/output files. To update FloPy from some branch of the [MODFLOW 6 repository](https://github.com/MODFLOW-USGS/modflow6), for instance the `develop` branch:

```shell
python -m flopy.mf6.utils.generate_classes --ref develop --no-backup
```

The `--repo` argument can be used to select an alternative repository (e.g. your fork). See [the FloPy documentation](https://flopy.readthedocs.io/en/stable/md/generate_classes.html) for more info.

## Running the examples

The example scripts can be run directly from the `scripts/` directory, e.g. `python ex-gwf-twri.py`. The environment variables described above can be used to control their behavior. The examples can also be run with `pytest`, converted to notebooks and/or executed with `jupytext`, or run as notebooks with `jupyter`.

**Note**: notebooks are *not* versioned in this repository &mdash; the `scripts/` are the single source of truth.

### Running with `pytest`

The examples can be tested from the `autotest/` directory, either directly as scripts, or as notebooks.

When run via `pytest`, behaviors can be controlled with `pytest` CLI flags:

- `--init`: just build models and write input files, defaults false
- `--no-write`: don't write models, defaults false
- `--no-run`: don't run models, defaults false
- `--plot`: enable plot creation, defaults false
- `--show`: show plots, defaults false
- `--no-save`: don't save static plots, defaults false
- `--no-gif`: don't create/save gifs, defaults false

The last three only apply if `--plot` is enabled. Plotting is disabled by default to avoid having to manually close plot widgets while running tests.

For instance, to create model input files (without running models or plotting results, to save time):

```shell
pytest -v -n auto test_scripts.py --init
```

### Running with `jupytext`

To convert an example script to a notebook manually with `jupytext` (add `--execute` to run the notebook after conversion):

```shell
jupytext --from py --to ipynb scripts/ex-gwf-twri.py -o notebooks/ex-gwf-twri.ipynb
```

### Running with `jupyter`

To start a Jupyter browser interface, run `jupyter notebook` from the `notebooks/` directory after notebooks have been created with `jupytext`.

## Contributing examples

To add a new example:

1. Add a new example script in the `scripts/` folder
2. Add a new LaTeX description in the `doc/sections/` folder
3. Ensure the script runs successfully and creates any expected files
3. Open a pull request from your fork to the upstream repo's `develop` branch

## Releasing the examples

GitHub Actions automatically creates a new release whenever code is merged into the `master` branch of this repository. Steps to prepare for a release include:

1. Run the examples to generate model input files, disabling model runs and plots &mdash; e.g. from the `autotest/` directory, run `pytest -v -n auto test_scripts.py`. **Note**: if model runs are enabled, the size of the `examples/` directory will balloon from double-digit MB to several GB. We only distribute input files, not output files.
2. Generate notebooks and tables &mdash; from the `scripts/` directory, run `process-scripts.py` to generate LaTeX tables for the documentation PDF from specially formatted code/comments in the example scripts.
3. Build the documentation PDF with multiple passes from e.g. `pdflatex` and `bibtex` &mdash; for instance, from the `doc/` directory:
    1. `pdflatex mf6examples.tex`
    2. `bibtex mf6examples`
    3. `pdflatex mf6examples.tex`
    4. `pdflatex mf6examples.tex`
4. Zip up the model input files.
5. Release the documentation PDF and model files archive.
