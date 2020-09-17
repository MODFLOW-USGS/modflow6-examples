import os
import re
import flopy

ex_pth = os.path.join("..", "examples")

def get_examples_list():
    # get order of examples from body.text
    ex_regex = re.compile("\\\\input{sections/(.*?)\\}")
    ex_order = []
    pth = os.path.join("..", "doc", "body.tex")
    with open(pth) as f:
        lines = f.read()
    for v in ex_regex.findall(lines):
        ex_order.append(v.replace(".tex", ""))

    # get list of all examples
    ex_list = []
    for name in sorted(os.listdir(ex_pth)):
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

def get_example_packages():
    ex_list = get_examples_list()
    ex_paks = {}
    for ex_name in ex_list:
        paks = []
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
                            if pak_type not in paks:
                                paks.append(pak_type)
                    for model_name in sim.model_names:
                        model = sim.get_model(model_name)
                        for pak in model.packagelist:
                            pak_type = pak.package_type
                            if pak_type not in paks:
                                paks.append(pak_type)
        # add packages for simulation to ex_paks dictionary
        ex_paks[ex_name] = paks
    return ex_paks

def build_tables(ex_paks):
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
                ex_list.append(ex_name)
        pak_dict[pak] = ex_list

    rtd_link = "https://modflow6.readthedocs.io/en/latest/_mf6io/"

    pth = os.path.join("..", ".doc", "introduction.md")
    f = open(pth, "w")

    line = "### Introduction\n\n"
    line += "This document describes example problems for MODFLOW 6. The examples " + \
            "demonstrate use of the capabilities for select components of " + \
            "MODFLOW 6. Examples have been included for the MODFLOW 6 components " + \
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

if __name__ == "__main__":
    ex_paks = get_example_packages()
    build_tables(ex_paks)