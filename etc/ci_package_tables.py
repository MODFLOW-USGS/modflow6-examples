import os
import sys
import re
import flopy

sys.path.append(os.path.join("..", "common"))
import build_table as bt

ex_pth = os.path.join("..", "examples")

def get_ordered_examples():
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

def get_examples_dict():
    ex_list = get_examples_list()
    ex_dict = {}
    for ex_name in ex_list:
        namefiles = []
        dimensions = []
        paks = []
        simulation_paks = []
        rootDir = os.path.join(ex_pth, ex_name)
        for dirName, subdirList, fileList in os.walk(rootDir):
            print('Found directory: {}'.format(dirName))
            for file_name in fileList:
                if file_name.lower() == "mfsim.nam":
                    print('\t{}'.format(file_name))
                    sim = flopy.mf6.MFSimulation.load(sim_ws=dirName, verbosity_level=0)
                    for pak in sim.sim_package_list:
                        if pak.package_abbr in ("gwfgwf", "gwfgwt",):
                            pak_type = pak.package_abbr
                            if pak_type not in simulation_paks:
                                simulation_paks.append(pak_type)
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
        paks = d["simulation_paks"]
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
                pak_link[ex_name] = "[{}](_examples/{}.html)".format(ex_name, ex_root)
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
    line += "This document describes example problems for MODFLOW 6. The examples " + \
            "demonstrate use of the capabilities for select components of " + \
            "MODFLOW 6. A pdf of the examples can be downloaded from " + \
            "[ReadtheDocs](https://modflow6-examples.readthedocs.io/en/latest/) " + \
            "or from the [current release](https://github.com/MODFLOW-USGS/" + \
            "modflow6-examples/releases/download/current/mf6examples.pdf) on " + \
            "[GitHub](https://github.com/MODFLOW-USGS/modflow6-examples/).\n\n" + \
            "Examples have been included for the MODFLOW 6 components " + \
            "summarized in the tables below.\n\n"
    f.write(line)

    join_fmt = ", "
    header = "| Package | Examples                                            " + \
             "               |\n" + "|---------|------------------------------" + \
             "--------------------------------------|\n"
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
    ex_order = get_ordered_examples()
    ex_tex = {}
    for idx, ex in enumerate(ex_order):
        for key, d in ex_dict.items():
            if ex in key:
                ex_number = [idx + 1] + [" " for i in range(len(d["paks"]) - 1)]
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
    col_widths = (0.10, 0.22, 0.25, 0.15, 0.28,)
    caption = "List of example problems and simulation characteristics."
    label = "tab:ex-table"

    lines = bt.get_header(caption, label, headings, col_widths=col_widths)

    for idx, (key, sim_dict) in enumerate(ex_tex.items()):
        for jdx, (ex_number, namefile, dimensions, paks) in enumerate(
                zip(
                    sim_dict["ex_number"],
                    sim_dict["namefiles"],
                    sim_dict["dimensions"],
                    sim_dict["paks"]
                    )
        ):
            if idx % 2 != 0:
                lines += "\t\t\\rowcolor{Gray}\n"
            lines += "\t\t"
            lines += "{} & ".format(ex_number)
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
    ex_dict = get_examples_dict()
    build_tex_tables(ex_dict)
    build_md_tables(ex_dict)
