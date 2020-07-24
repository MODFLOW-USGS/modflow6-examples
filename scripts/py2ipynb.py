import os

tpth = "raw.py"

# only process python files starting with ex_
files = [
    file
    for file in os.listdir()
    if file.endswith(".py") and file.startswith("ex_")
]


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
            modifyIndent = len(lines[idx+1]) - len(lines[idx+1].lstrip(' '))
            continue

        # exclude nosetest functions
        if "# nosetest" in line:
            if skip:
                skip = False
                continue
            else:
                skip = True
        if skip:
            continue
        f.write("{}\n".format(line[modifyIndent:]))
    f.close()

    #
    basename = os.path.splitext(file)[0] + ".ipynb"
    opth = os.path.join("..", "notebooks", basename)
    os.system("p2j {} -o -t {}".format(tpth, opth))

    # remove temporary file
    if os.path.isfile(tpth):
        os.remove(tpth)
