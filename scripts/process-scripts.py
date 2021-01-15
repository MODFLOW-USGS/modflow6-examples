import os
import sys
import re
import flopy

sys.path.append(os.path.join("..", "common"))
import build_table as bt

# path to the example files
ex_pth = os.path.join("..", "examples")

# only process python files starting with ex_
files = [
    file
    for file in os.listdir()
    if file.endswith(".py") and file.startswith("ex-")
]


def _replace_quotes(proc_str):
    for r in (("'", ""), ('"', "")):
        proc_str = proc_str.replace(*r)
    return proc_str


def make_notebooks():
    nb_pth = os.path.join("..", "notebooks")
    if not os.path.isdir(nb_pth):
        os.makedirs(nb_pth)

    tpth = "raw.py"
    # converts python scripts to jupyter notebooks
    for file in files:
        # do a little processing
        with open(file) as f:
            lines = f.read().splitlines()
        f = open(tpth, "w")
        skip = False
        modifyIndent = 0
        for idx, line in enumerate(lines):
            # exclude if __name__ == "main"
            if "if __name__" in line:
                modifyIndent = len(lines[idx + 1]) - len(
                    lines[idx + 1].lstrip(" ")
                )
                continue

            # exclude nosetest functions
            if "# nosetest" in line.lower():
                if skip:
                    skip = False
                    continue
                else:
                    skip = True
            if skip:
                continue
            f.write("{}\n".format(line[modifyIndent:]))
        f.close()

        # convert temporary python file to a notebook
        basename = os.path.splitext(file)[0] + ".ipynb"
        opth = os.path.join(nb_pth, basename)
        cmd = (
            "jupytext",
            "--from py",
            "--to ipynb",
            "-o",
            opth,
            tpth,
        )
        print(" ".join(cmd))
        os.system(" ".join(cmd))

        # remove temporary file
        if os.path.isfile(tpth):
            os.remove(tpth)


def table_standard_header(caption, label):
    col_widths = (
        0.5,
        0.3,
    )
    headings = (
        "Parameter",
        "Value",
    )
    return bt.get_header(
        caption, label, headings, col_widths=col_widths, center=False
    )


def table_scenario_header(caption, label):
    col_widths = (
        0.1,
        0.25,
        0.3,
        0.15,
    )
    headings = (
        "Scenario",
        "Scenario Name",
        "Parameter",
        "Value",
    )
    return bt.get_header(
        caption, label, headings, col_widths=col_widths, center=False
    )


def table_footer():
    return bt.get_footer()


def make_tables():
    tab_pth = os.path.join("..", "tables")
    if not os.path.isdir(tab_pth):
        os.makedirs(tab_pth)

    for file in files:
        print("processing...'{}'".format(file))
        basename = os.path.splitext(os.path.basename(file))[0].replace(
            "_", "-"
        )
        # do a little processing
        with open(file) as f:
            lines = f.read().splitlines()

        # parse optional parameters dictionary
        parameters = None
        tag = "parameters ="
        for idx, line in enumerate(lines):
            if line.lower().startswith(tag):
                dict_str = line[len(tag) + 1 :].strip()
                for jdx in range(idx + 1, len(lines)):
                    if len(lines[jdx].strip()) < 1:
                        break
                    dict_str += " " + lines[jdx].strip()

                # parse the dictionary string into a dictionary
                parameters = eval(dict_str)
                break

        # parse optional parameter units dictionary, if parameters are
        # specified in the script
        if parameters is not None:
            parameter_units = None
            tag = "parameter_units ="
            for idx, line in enumerate(lines):
                if line.lower().startswith(tag):
                    dict_str = line[len(tag) + 1 :].strip()
                    for jdx in range(idx + 1, len(lines)):
                        if len(lines[jdx].strip()) < 1:
                            break
                        dict_str += " " + lines[jdx].strip()

                    # parse the dictionary string into a dictionary
                    parameter_units = eval(dict_str)
                    break

        # create scenario table if parameters are specified in the script
        if parameters is not None:
            tab_name = "{}-scenario".format(basename)
            caption = "Scenario parameters for example {}.".format(basename)
            label = "tab:{}".format(tab_name)
            fpth = os.path.join(tab_pth, tab_name + ".tex")
            f = open(fpth, "w")

            f.write(table_scenario_header(caption, label))

            scenario_count = 0
            for scenario_name, value_dict in parameters.items():
                if scenario_count % 2 != 0:
                    row_color = "\t\t\\rowcolor{Gray}\n"
                else:
                    row_color = ""
                scenario_count += 1
                table_line = "{} & {} & ".format(scenario_count, scenario_name)
                for text, value in value_dict.items():
                    text_to_write = text.replace("_", "~")
                    units = ""
                    if parameter_units is not None:
                        try:
                            units = " ({})".format(parameter_units[text])
                        except:
                            units = " (unknown)"
                    if len(table_line) > 0:
                        table_line += "\t{}{} & {}".format(
                            text_to_write, units, value
                        )
                    else:
                        table_line = "\t& & {}{} & {}".format(
                            text_to_write, units, value
                        )
                    table_line += " \\\\\n"
                    if len(row_color) > 0:
                        f.write(row_color)
                    f.write(table_line)
                    table_line = ""

            # write footer
            f.write(table_footer())

            # finalize table
            f.close()

        # tables
        table_number = 0
        for idx, line in enumerate(lines):
            table_text = []
            table_value = []
            if line.lower().startswith("# table"):
                for table_line in lines[idx + 1 :]:
                    # skip empty lines
                    if len(table_line.strip()) < 1:
                        continue
                    if table_line.startswith("# "):
                        break
                    ipos = table_line.find("# ")
                    if ipos > 0:
                        table_text.append(table_line[ipos + 1 :].strip())
                        table_value.append(
                            _replace_quotes(
                                table_line[0:ipos].split("=")[1].strip()
                            )
                        )
            if len(table_text) > 0:
                table_number += 1
                tab_name = "{}-{:02d}".format(basename, table_number)
                caption = "Model parameters for example {}.".format(basename)
                label = "tab:{}".format(tab_name)
                fpth = os.path.join(tab_pth, tab_name + ".tex")
                f = open(fpth, "w")

                f.write(table_standard_header(caption, label))

                # write table
                line_count = 0
                for text, value in zip(table_text, table_value):
                    if line_count % 2 != 0:
                        f.write("\t\\rowcolor{Gray}\n")
                    f.write("\t{} & {} \\\\\n".format(text, value))
                    line_count += 1

                # write footer
                f.write(table_footer())

                # finalize table
                f.close()


def get_ordered_examples():
    print("creating a ordered list of examples from the LaTeX document")
    # get order of examples from body.text
    ex_regex = re.compile("\\\\input{sections/(.*?)\\}")
    ex_order = []
    pth = os.path.join("..", "doc", "body.tex")
    with open(pth) as f:
        lines = f.read()
    for v in ex_regex.findall(lines):
        ex_order.append(v.replace(".tex", ""))
    return ex_order


def get_examples_list():
    print("creating a list of available examples")
    # examples to exclude
    exclude = ("ex-gwf-csub-p02c",)

    # get order of examples from body.text
    ex_order = get_ordered_examples()

    # get list of all examples
    ex_list = []
    for name in sorted(os.listdir(ex_pth)):
        if name in exclude:
            continue
        pth = os.path.join(ex_pth, name)
        if os.path.isdir(pth):
            ex_list.append(name)

    # create final list with all of the examples in body.tex order
    ex_final = []
    for o in ex_order:
        for n in ex_list:
            if n.startswith(o):
                if n not in ex_final:
                    ex_final.append(n)

    # add any orphans
    if len(ex_final) != len(ex_list):
        for n in ex_list:
            if n not in ex_final:
                ex_final.append(n)

    return ex_final


def get_examples_dict(verbose=False):
    print("creating dictionary with examples information")
    ex_list = get_examples_list()
    ex_dict = {}
    for ex_name in ex_list:
        namefiles = []
        dimensions = []
        paks = []
        simulation_paks = []
        rootDir = os.path.join(ex_pth, ex_name)
        for dirName, subdirList, fileList in os.walk(rootDir):
            if verbose:
                print("Found directory: {}".format(dirName))
            for file_name in fileList:
                if file_name.lower() == "mfsim.nam":
                    if verbose:
                        msg = (
                            "  Found MODFLOW 6 simulation "
                            + "name file: {}".format(file_name)
                        )
                        print(msg)
                    sim = flopy.mf6.MFSimulation.load(
                        sim_ws=dirName, verbosity_level=0
                    )
                    sim_paks = []
                    for pak in sim.sim_package_list:
                        pak_type = pak.package_abbr
                        if pak_type not in simulation_paks:
                            sim_paks.append(pak_type)
                    for model_name in sim.model_names:
                        model = sim.get_model(model_name)
                        namefiles.append(model.namefile)
                        dimensions.append(model.modelgrid.shape)
                        model_paks = []
                        for pak in model.packagelist:
                            pak_type = pak.package_type
                            if pak_type not in model_paks:
                                model_paks.append(pak_type)
                                if pak_type == "npf":
                                    if model.npf.xt3doptions.array:
                                        model_paks.append("xt3d")
                        simulation_paks.append(sim_paks)
                        paks.append(model_paks)
        # add packages for simulation to ex_paks dictionary
        ex_dict[ex_name] = {
            "namefiles": namefiles,
            "dimensions": dimensions,
            "paks": paks,
            "simulation_paks": simulation_paks,
        }
    return ex_dict


def build_md_tables(ex_dict):
    print("building markdown tables for ReadtheDocs")
    # determine the order of the examples from the LaTeX document
    ex_order = get_ordered_examples()

    # create a dictionary in the correct order
    ex_md = {}
    for ex in ex_order:
        for key, d in ex_dict.items():
            if ex in key:
                ex_md[key] = d

    # build a dictionary with all of the unique packages for a example
    ex_paks = {}
    for key, d in ex_md.items():
        paks = d["simulation_paks"][0]
        # add model packages that have not ben added yet
        for model_paks in d["paks"]:
            # this is the union of the two lists
            paks = list(set(paks) | set(model_paks))
        ex_paks[key] = paks

    # build dictionary with hyperlinks
    pak_link = {}
    for ex_name in ex_paks.keys():
        for ex_root in ex_order:
            if ex_root in ex_name:
                pak_link[ex_name] = "[{}](_examples/{}.html)".format(
                    ex_name, ex_root
                )
                break
        if ex_name not in list(pak_link.keys()):
            pak_link[ex_name] = ex_name

    # build list of unique packages
    pak_unique = []
    for ex_name, paks in ex_paks.items():
        for pak in paks:
            if pak not in pak_unique:
                pak_unique.append(pak)

    pak_dict = {}
    for pak in pak_unique:
        ex_list = []
        for ex_name, paks in ex_paks.items():
            if pak in paks:
                ex_list.append(pak_link[ex_name])
        pak_dict[pak] = ex_list

    rtd_link = "https://modflow6.readthedocs.io/en/latest/_mf6io/"

    pth = os.path.join("..", ".doc", "introduction.md")
    f = open(pth, "w")

    line = "### Introduction\n\n"
    line += (
        "This document describes example problems for MODFLOW 6. The examples "
        + "demonstrate use of the capabilities for select components of "
        + "MODFLOW 6. A pdf of the examples can be downloaded from "
        + "[ReadtheDocs](https://modflow6-examples.readthedocs.io/en/latest/) "
        + "or from the [current release](https://github.com/MODFLOW-USGS/"
        + "modflow6-examples/releases/download/current/mf6examples.pdf) on "
        + "[GitHub](https://github.com/MODFLOW-USGS/modflow6-examples/).\n\n"
        + "Examples have been included for the MODFLOW 6 components "
        + "summarized in the tables below.\n\n"
    )
    f.write(line)

    join_fmt = ", "
    header = (
        "| Package | Examples                                            "
        + "               |\n"
        + "|---------|------------------------------"
        + "--------------------------------------|\n"
    )
    footer = "\n\n"

    pak_link = {
        "dis": "gwf-dis.html",
        "disv": "gwf-disv.html",
        "disu": "gwf-disu.html",
    }
    line = "#### Discretization types\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    pak_link = {
        "gwfgwf": "exg-gwfgwf.html",
        "gwfgwt": "exg-gwfgwt.html",
    }
    line = "#### Exchanges\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    pak_link = {
        "npf": "gwf-npf.html",
        "hfb": "gwf-hfb.html",
        "sto": "gwf-sto.html",
        "csub": "gwf-csub.html",
    }
    line = "#### Groundwater Flow Model Internal Flow Packages\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    pak_link = {
        "rch": "gwf-rch.html",
        "rcha": "gwf-rcha.html",
        "evt": "gwf-evt.html",
        "evta": "gwf-evta.html",
        "chd": "gwf-chd.html",
        "drn": "gwf-drn.html",
        "ghb": "gwf-ghb.html",
        "riv": "gwf-riv.html",
        "wel": "gwf-wel.html",
        "buy": "gwf-buy.html",
    }
    line = "#### Groundwater Flow Model Standard Boundary Packages\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    pak_link = {
        "lak": "gwf-lak.html",
        "sfr": "gwf-sfr.html",
        "maw": "gwf-maw.html",
        "uzf": "gwf-uzf.html",
        "mvr": "gwf-mvr.html",
    }
    line = "#### Groundwater Flow Model Advanced Boundary Packages\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    pak_link = {
        "adv": "gwt-adv.html",
        "dsp": "gwt-dsp.html",
        "mst": "gwt-mst.html",
        "ist": "gwt-ist.html",
        "fmi": "gwt-fmi.html",
    }
    line = "#### Groundwater Transport Model Internal Flow Packages\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer + "\n\n"
    f.write(line)

    pak_link = {
        "ssm": "gwt-ssm.html",
        "src": "gwt-src.html",
        "cnc": "gwt-cnc.html",
    }
    line = "#### Groundwater Transport Model Standard Boundary Packages\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    pak_link = {
        "lkt": "gwt-lkt.html",
        "sft": "gwt-sft.html",
        "mwt": "gwt-mwt.html",
        "uzt": "gwt-uzt.html",
        "mvt": "gwt-mvt.html",
    }
    line = "#### Groundwater Transport Model Advanced Boundary Packages\n\n"
    line += header
    for pak, value in pak_link.items():
        if pak in pak_dict.keys():
            line += "| [{}]({}{}) |".format(pak.upper(), rtd_link, value)
            line += " {} |\n".format(join_fmt.join(pak_dict[pak]))
    line += footer
    f.write(line)

    f.close()


def build_tex_tables(ex_dict):
    print("building LaTeX table for mf6examples.tex introduction")
    ex_order = get_ordered_examples()
    ex_tex = {}
    for idx, ex in enumerate(ex_order):
        for key, d in ex_dict.items():
            if ex in key:
                ex_number = [idx + 1] + [
                    " " for i in range(len(d["paks"]) - 1)
                ]
                d["ex_number"] = ex_number
                ex_tex[key] = d

    # build latex table for pdf document
    headings = (
        "Example",
        "Simulation",
        "Namefile(s)",
        "Grid \\newline Dimensions",
        "Packages",
    )
    col_widths = (
        0.10,
        0.22,
        0.25,
        0.15,
        0.28,
    )
    caption = "List of example problems and simulation characteristics."
    label = "tab:ex-table"

    lines = bt.get_header(
        caption, label, headings, col_widths=col_widths, firsthead=True
    )

    on_ex = 0
    for idx, (key, sim_dict) in enumerate(ex_tex.items()):
        for jdx, (
            ex_number,
            namefile,
            dimensions,
            model_paks,
            sim_paks,
        ) in enumerate(
            zip(
                sim_dict["ex_number"],
                sim_dict["namefiles"],
                sim_dict["dimensions"],
                sim_dict["paks"],
                sim_dict["simulation_paks"],
            )
        ):
            # union of sim_paks and model_paks to create unique list of packages
            paks = sim_paks
            for pak in model_paks:
                if pak not in paks:
                    paks.append(pak)

            # eliminate any duplicate packages
            paks = sorted(list(set(paks)))

            # separate examples by a line
            if isinstance(ex_number, int):
                if ex_number > on_ex:
                    on_ex = ex_number
                    include_idx = True
                    lines += "\t\t\\hline\n"
            if on_ex % 2 == 0:
                lines += "\t\t\\rowcolor{Gray}\n"
            lines += "\t\t"
            if include_idx:
                include_idx = False
                lines += "{} & ".format(ex_number)
            else:
                lines += " & "
            if jdx == 0:
                lines += "{} & ".format(key)
            else:
                lines += " & "
            lines += " {} & ".format(namefile.replace("_", "\\_"))
            lines += " {} & ".format(dimensions)

            # packages
            pak_line = []
            for pak in paks:
                pak_line.append(pak.upper())
            lines += " ".join(pak_line) + " \\\\\n"

    lines += bt.get_footer()

    # create table
    pth = os.path.join("..", "tables", "ex-table.tex")
    f = open(pth, "w")
    f.write(lines)
    f.close()


if __name__ == "__main__":
    verbose = False
    # parse command line arguments
    for arg in sys.argv:
        if arg in ("-v", "--verbose"):
            verbose = True

    # make the notebooks
    make_notebooks()

    # make the tables from the scripts
    make_tables()

    # get a dictionary with example information from the example files
    ex_dict = get_examples_dict(verbose=verbose)

    # make the summary table for the LaTeX document
    build_tex_tables(ex_dict)

    # make the markdown table for ReadtheDocs
    build_md_tables(ex_dict)
