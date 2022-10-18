# Creating New MODFLOW 6 Examples

Adding a new example requires adding a new example script in the `scripts` folder and adding a new latex problem description in the `doc/sections` folder.  The following steps describe the process for running all the scripts, converting them to jupyter notebooks and latex tables, and then building the pdf document.

## Installation of Modeling Programs

Install modeling programs into the `bin` directory.  This can be done using the get-modflow utility that is available when the flopy Python package is installed.  From the `bin` directory, issue the following commands:

```commandline
get-modflow .
get-modflow --repo modflow6-nightly-build .
```

There are quite a few Python packages that are required in order to create MODFLOW 6 examples.  They are listed in `etc/requirements.pip.txt` and `etc/requirements.usgs.txt`.

## Update FloPy

It's important that FloPy is up-to-date with the latest changes to MODFLOW 6.  This can be achieved by running the following Python commands:

```commandline
import flopy
flopy.mf6.utils.generate_classes(branch="develop", backup=False)
```

## Running the Example Scripts

Pytest can be used to run all of the examples.  From the `scripts` directory, issue the following command:

```commandline
pytest -v -n=auto ex-*.py 
```

or use the following command, which runs in serial instead of running in parallel:

```commandline
pytest -v ex-*.py 
```

## Extract Notebooks and Tables from Scripts

The example scripts are converted into jupyter notebooks, and latex tables are created from the scripts using jupytext.  To convert all of the scripts, run the following command in the scrips directory:

```commandline
python process-scripts.py
```

## Build the Latex Document

If the figures and tables were generated correctly, then it should be possible to build the pdf document for the examples.  The pdf is created by converting the `doc/mf6examples.tex` document into a pdf.