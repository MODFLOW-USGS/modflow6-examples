import os

# only process python files starting with ex_
files = [
    file for file in os.listdir() if file.endswith(".py") and file.startswith("ex-")
]


def _replace_quotes(proc_str):
    for r in (("'", ""), ('"', "")):
        proc_str = proc_str.replace(*r)
    return proc_str


def make_notebooks():
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
                modifyIndent = len(lines[idx + 1]) - len(lines[idx + 1].lstrip(" "))
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
        opth = os.path.join("..", "notebooks", basename)
        os.system("p2j {} -o -t {}".format(tpth, opth))

        # remove temporary file
        if os.path.isfile(tpth):
            os.remove(tpth)


def make_tables():
    for file in files:
        basename = os.path.splitext(os.path.basename(file))[0].replace("_", "-")
        # do a little processing
        with open(file) as f:
            lines = f.read().splitlines()

        # parameters table
        for idx, line in enumerate(lines):
            if line.lower().startswith("parameters ="):
                dict_str = line[12:].strip()
                for jdx in range(idx + 1, len(lines)):
                    if len(lines[jdx].strip()) < 1:
                        break
                    dict_str += " " + lines[jdx].strip()

                # parse the dictionary string into a dictionary
                parameters = eval(dict_str)

                # scenario table
                tab_name = "{}-scenario".format(basename)
                fpth = os.path.join("..", "tables", tab_name + ".tex")
                f = open(fpth, "w")

                scenario_count = 0
                for scenario_name, value_dict in parameters.items():
                    if scenario_count % 2 != 0:
                        row_color = "\\rowcolor{Gray}\n"
                    else:
                        row_color = ""
                    scenario_count += 1
                    table_line = "{} & {} & ".format(scenario_count, scenario_name)
                    for text, value in value_dict.items():
                        if len(table_line) > 0:
                            table_line += "{} & {}".format(text, value)
                        else:
                            table_line = "& & {} & {}".format(text, value)
                        table_line += " \\\\\n"
                        if len(row_color) > 0:
                            f.write(row_color)
                        f.write(table_line)
                        table_line = ""

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
                            _replace_quotes(table_line[0:ipos].split("=")[1].strip())
                        )
            if len(table_text) > 0:
                table_number += 1
                tab_name = "{}{:02d}".format(basename, table_number)
                fpth = os.path.join("..", "tables", tab_name + ".tex")
                f = open(fpth, "w")

                # write table
                line_count = 0
                for text, value in zip(table_text, table_value):
                    if line_count % 2 != 0:
                        f.write("\\rowcolor{Gray}\n")
                    f.write("{} & {} \\\\\n".format(text, value))
                    line_count += 1

                # finalize table
                f.close()


if __name__ == "__main__":
    make_notebooks()
    make_tables()
