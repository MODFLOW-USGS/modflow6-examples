<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [MODFLOW 6 Example Problems](#modflow-6-example-problems)
  - [Release contents](#release-contents)
  - [Repository contents](#repository-contents)
  - [Resources](#resources)
  - [Issues](#issues)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# MODFLOW 6 Example Problems

[![CI](https://github.com/MODFLOW-USGS/modflow6-examples/actions/workflows/ex-workflow.yml/badge.svg)](https://github.com/MODFLOW-USGS/modflow6-examples/actions/workflows/ex-workflow.yml)
[![Documentation Status](https://readthedocs.org/projects/modflow6-examples/badge/?version=latest)](https://modflow6-examples.readthedocs.io/en/latest/?badge=latest)

This repository contains MODFLOW 6 example problems.

## Release contents

When changes reach the `master` branch, this repository's contents are [rebuilt and posted as a new release](https://github.com/MODFLOW-USGS/modflow6-examples/releases) with the following assets:

* PDF containing a description of the example problems
* archive containing input files for the example models

## Repository contents

* `autotest/`: tests for the example models
* `data/`: data files used by the example models
* `doc/`: LaTeX files to build the PDF documentation
* `etc/`: dependency specifications and miscellaneous scripts
* `images/`: static images used in the documentation
* `scripts/`: Python scripts that use FloPy to create, run, and post-process the example models
* `.github/`: continuous integration workflows
* `.doc/`: source files for the ReadTheDocs site

## Resources

* [U.S. Geological Survey MODFLOW 6 Page](https://www.usgs.gov/software/modflow-6-usgs-modular-hydrologic-model)
* [MODFLOW 6 binary executables for Windows, Mac, and Linux](https://github.com/MODFLOW-USGS/modflow6-nightly-build/releases)
* [MODFLOW 6 Repository](https://github.com/MODFLOW-USGS/modflow6)
* [Binary executables for MODFLOW and other related programs that may be required to run these examples](https://github.com/MODFLOW-USGS/executables)

## Issues

Any issues with the MODFLOW 6 example problems should be posted on the main [MODFLOW 6 GitHub repo](https://github.com/MODFLOW-USGS/modflow6) and flagged with the [examples](https://github.com/MODFLOW-USGS/modflow6/labels/examples) label.
