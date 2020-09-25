import os
import sys

sys.path.append(os.path.join("..", "common"))
import build_table as bt

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
    col_widths = (0.5, 0.3,)
    headings = (
        "Parameter",
        "Value",
    )
    return bt.get_header(caption, label, headings,
                         col_widths=col_widths, center=False)

def table_scenario_header(caption, label):
    col_widths = (0.1, 0.25, 0.3, 0.15, )
    headings = (
        "Scenario",
        "Scenario Name",
        "Parameter",
        "Value",
    )
    return bt.get_header(caption, label, headings,
                         col_widths=col_widths, center=False)

def table_footer():
    return bt.get_footer()

def make_tables():
    tab_pth = os.path.join("..", "tables")
    if not os.path.isdir(tab_pth):
        os.makedirs(tab_pth)

    for file in files:
        print("processing...'{}'".format(file))
        basename = os.path.splitext(os.path.basename(file))[0].replace("_", "-")
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
            label = "tab:{}".format(tab_name)
            fpth = os.path.join(tab_pth, tab_name + ".tex")
            f = open(fpth, "w")

            f.write(table_scenario_header("Model Scenario Parameters", label))

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
                label = "tab:{}".format(tab_name)
                fpth = os.path.join(tab_pth, tab_name + ".tex")
                f = open(fpth, "w")

                f.write(table_standard_header("Model Parameters", label))

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


if __name__ == "__main__":
    make_notebooks()
    make_tables()
