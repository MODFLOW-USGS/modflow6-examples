import os
import re
import shutil
from subprocess import Popen, PIPE

# check that pandoc is available and the version is sufficient
args = ("pandoc", "--version")
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=".")
stdout, stderr = proc.communicate()
msg = ""
if stdout:
    for line in stdout.decode("utf-8").splitlines():
        if "pandoc" in line:
            t = line.split()[1]
            version = []
            for v in t.split("."):
                version.append(int(v))
            version = tuple(version)
            if version[0:2] < (2, 11):
                msg = "pandoc 2.11 or higher must be available"
            break
if stderr:
    msg = stderr.decode("utf-8")
if len(msg) > 0:
    raise RuntimeError(msg)
else:
    print("\npandoc version: {}\n".format(t))

# get list of examples from body.text
ex_regex = re.compile("\\\\input{sections/(.*?)\\}")
ex_list = []
pth = os.path.join("..", "doc", "body.tex")
with open(pth) as f:
    lines = f.read()
for v in ex_regex.findall(lines):
    ex_list.append(v.replace(".tex", ""))

# get list of example problems
ex_dirs = []
ex_pth = os.path.join("..", "examples")
for name in sorted(os.listdir(ex_pth)):
    pth = os.path.join(ex_pth, name)
    if os.path.isdir(pth):
        ex_dirs.append(name)

# create examples rst
pth = os.path.join("..", ".doc", "examples.rst")
print("creating...'{}'".format(pth))
f = open(pth, "w")
line = "============================\n"
line += "MODFLOW 6 – Example problems\n"
line += "============================\n\n\n"
line += ".. toctree::\n"
line += "   :numbered:\n"
line += "   :maxdepth: 1\n\n"
for ex in ex_list:
    line += "   _examples/{}.rst\n".format(ex)
f.write(line)
f.close()

# create rtd examples directory
dst = os.path.join("..", ".doc", "_examples")
print("cleaning and creating...'{}'".format(dst))
if os.path.isdir(dst):
    shutil.rmtree(dst)
os.makedirs(dst, exist_ok=True)

# read base latex file
pth = os.path.join("..", "doc", "mf6examples.tex")
print("reading...'{}'".format(pth))
with open(pth) as f:
    orig_latex = f.readlines()

latex_tag = "\\input{./body.tex}"
frontmatter_tag = "\\input{./frontmatter.tex}"
doc_pth = os.path.join("..", "doc")
for ex in ex_list:
    print("creating restructured text file for {} example".format(ex))
    src = os.path.join("..", "doc", "ex.tex")
    f = open(src, "w")
    for line in orig_latex:
        if latex_tag in line:
            new_tag = "\\input{{sections/{}.tex}}".format(ex)
            line = line.replace(latex_tag, new_tag)
        elif frontmatter_tag in line:
            line = line.replace(frontmatter_tag, "")
        f.write(line)
    f.close()

    # create restructured text file for example using using pandoc
    dst = os.path.join("..", ".doc", "_examples", "{}.rst".format(ex))
    print("running pandoc to create {}".format(dst))
    args = (
        "pandoc",
        "-s",
        "--citeproc",
        "-f",
        "latex",
        "-t",
        "rst",
        "--bibliography=mf6examples.bib",
        "--csl=wrr.csl",
        "ex.tex",
        "-o",
        dst,
    )

    print(" ".join(args))
    proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=doc_pth)
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print("Errors:\n{}".format(stderr.decode("utf-8")))

    # read restructured text file as a string
    print("reading...'{}'".format(dst))
    with open(dst, "r") as file:
        lines = file.read()

    # find equation labels in lines
    ex_regex = re.compile("\\\\label{(.*?)\\}")
    replace_eq_labels = {}
    eq_labels = []
    for v in ex_regex.findall(lines):
        tag = "\\label{{{}}}".format(v)
        label = "   :label: {}".format(v.replace(":", "-"))
        eq_labels.append(label)
        tag = "`[{0}] <#{0}>`__".format(v)
        replace_eq_labels[tag] = ":eq:`{}`".format(v.replace(":", "-"))

    # find figure references in lines
    ex_regex = re.compile("\\`(.*?)\\`__")
    ex_tag = re.compile("\\#(.*?)\\>")
    fig_tab_refs = {}
    for v in ex_regex.findall(lines):
        tag0 = "`{}`__".format(v)
        for tag in ex_tag.findall(tag0):
            if tag0 not in list(fig_tab_refs.keys()):
                fig_tab_refs[tag0] = ":numref:`{}`".format(tag)

    # read restructured text file for example
    print("reading...'{}'".format(dst))
    with open(dst) as f:
        lines = f.readlines()

    # editing restructured text file for example
    print("editing...'{}'".format(dst))
    f = open(dst, "w")

    write_line = True
    in_reference = False
    eq_no = 0
    in_table = False
    table_name = None
    for idx, line in enumerate(lines):
        # skip the title
        if idx < 6:
            continue

        # replace latex equation labels
        for key, value in replace_eq_labels.items():
            if key in line:
                line = line.replace(key, value)

        # replace figure and table references
        for key, value in fig_tab_refs.items():
            if key in line:
                line = line.replace(key, value)

        tag = " "
        if tag in line:
            line = line.replace(tag, " ")

        tag = "tab:ex-"
        if tag in line:
            line = line.replace(tag, "tab-ex-")

        tag = "fig:ex-"
        if tag in line:
            line = line.replace(tag, "fig-ex-")

        tag = "\\label{"
        if tag in line:
            continue

        tag = ".. math::"
        if tag in line:
            line = "{}\n{}\n".format(tag, eq_labels[eq_no])
            eq_no += 1

        tag = ".. figure:: ../figures/"
        if tag in line:
            line = line.replace(tag, ".. figure:: ../_images/")

        tag = ".. figure:: ../images/"
        if tag in line:
            line = line.replace(tag, ".. figure:: ../_images/")

        tag = ":alt:"
        if tag in line:
            write_line = False

        tag = ":name:"
        if not write_line and tag in line:
            write_line = True

        # table tabs
        tags = (
            "Table :numref:",
            "table :numref:",
            "Tables :numref:",
        )
        for tag in tags:
            if tag in line:
                line = line.replace(tag, ":numref:")

        # figure tags
        tags = (
            "fig. :numref:",
            "Figure :numref:",
            "figure :numref:",
            "figures :numref:",
        )
        for tag in tags:
            if tag in line:
                line = line.replace(tag, ":numref:")

        tag = ".. container::"
        if tag == line.strip():
            in_table = True
            continue

        if in_table:
            if len(line.lstrip()) == len(line):
                line = "\n" + line
                in_table = False

        tag = ":name: tab-ex-"
        if tag in line:
            table_name = line.split()[1]
            continue

        tag = ".. table::"
        if tag in line:
            line = ".. _{}:\n\n".format(table_name) + line.lstrip() + "\n"

        tag = ".. container:: references"
        if tag in line:
            in_reference = True
            line = "References Cited\n----------------\n\n"

        if in_reference:
            tag = ".. container::"
            if tag in line:
                continue

            tag = ":name:"
            if tag in line:
                continue

            if line.startswith("      "):
                line = line.lstrip()

        if write_line:
            f.write(line)

    # Jupyter Notebook section
    ex_root = ex.replace(".tex", "")

    # write Jupyter Notebook section
    line = "\n\nJupyter Notebook\n"
    line += "----------------\n\n"
    line += (
        "The Jupyter notebook used to create the MODFLOW 6 input files\n"
        + "for this example and post-process the results is:\n\n"
    )
    line += "* `{0} <../_notebooks/{0}.html>`_\n".format(ex_root)
    line += "\n"

    # Check to see if there is a gif with the same name as the example name
    # (e.g. ex-gwt-saltlake.gif) and if so, then assume this is an animated
    # gif and add an Animation section to the restructured text file
    fname = "{}.gif".format(ex)
    if fname in os.listdir("../figures"):
        line = "\n\nAnimation\n"
        line += "----------------\n\n"
        line += "Animation of model results:\n\n"
        line += ".. image:: ../_images/{}".format(fname)
    line += "\n"

    # Write the restructured text lines and close the file
    f.write(line)
    f.close()

    # # convert the rst tables
    # print("running pandoc to create {}".format(dst))
    # args = (
    #     "rstlisttable",
    #     dst,
    # )
    # print(" ".join(args))
    # proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=doc_pth)
    # stdout, stderr = proc.communicate()
    # if stderr:
    #     print("Errors:\n{}".format(stderr.decode("utf-8")))
    # else:
    #     f = open(dst, "w")
    #     f.write(stdout.decode("utf-8"))
    #     f.close()

    # clean up temporary latex file
    if os.path.isfile(src):
        os.remove(src)

# create rtd figure directory
dst = os.path.join("..", ".doc", "_images")
print("cleaning and creating...'{}'".format(dst))
if os.path.isdir(dst):
    shutil.rmtree(dst)
os.makedirs(dst, exist_ok=True)

# copy figures to rtd directory
src_dirs = (os.path.join("..", "figures"), os.path.join("..", "images"))
for src_dir in src_dirs:
    file_names = [
        file_name
        for file_name in os.listdir(src_dir)
        if os.path.isfile(os.path.join(src_dir, file_name))
        and (file_name.endswith(".png") or file_name.endswith(".gif"))
    ]
    for file_name in file_names:
        src = os.path.join(src_dir, file_name)
        print("copy '{}' -> '{}' directory".format(src, dst))
        shutil.copy2(src, dst)
