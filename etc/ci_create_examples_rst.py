import os
import shutil
from subprocess import Popen, PIPE

# create examples.rst using pandoc
print("running pandoc to create examples.rst")
doc_pth = os.path.join("..", "doc")
args = (
    "pandoc",
    "-s",
    "-f",
    "latex",
    "-t",
    "rst",
    "--bibliography=mf6examples.bib",
    "mf6examples.tex",
    "-o",
    "../.doc/examples.rst"
)
print(" ".join(args))
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=doc_pth)
stdout, stderr = proc.communicate()
if stdout:
    print(stdout.decode("utf-8"))
if stderr:
    print("Errors:\n{}".format(stderr.decode("utf-8")))

# read examples.rst
pth = os.path.join("..", ".doc", "examples.rst")
print("reading...'{}'".format(pth))
with open(pth) as f:
    lines = f.readlines()

# editing examples.rst
print("editing...'{}'".format(pth))
f = open(pth, "w")

# find section headings
heading_pos = []
line_len0 = len(lines[0].strip())
for idx, line in enumerate(lines):
    # write the title
    if idx < 4:
        f.write(line)
        continue

    # determine if this is a section heading
    tag = line_len0 * "="
    if line.strip() == tag:
        heading_pos.append(idx - 1)

    line_len0 = len(line.strip())

istart = heading_pos[1]
write_line = True
for line in lines[istart:]:
    tag = ".. figure:: ../figures/"
    if tag in line:
        line = line.replace(tag, ".. figure:: _images/")

    tag = ".. figure:: ../images/"
    if tag in line:
        line = line.replace(tag, ".. figure:: _images/")

    tag = ":alt:"
    if tag in line:
        write_line = False

    tag = ":name:"
    if not write_line and tag in line:
        write_line = True

    if write_line:
        f.write(line)

f.close()

# create rtd figure directory
dst = os.path.join("..", ".doc", "_images")
print("cleaning and creating...'{}'".format(dst))
if os.path.isdir(dst):
    shutil.rmtree(dst)
os.makedirs(dst)

# copy figures to rtd directory
src_dirs = (os.path.join("..", "figures"),
            os.path.join("..", "images"))
for src_dir in src_dirs:
    file_names = [file_name for file_name in os.listdir(src_dir) if os.path.isfile(
        os.path.join(src_dir, file_name)) and file_name.endswith(".png")]
    for file_name in file_names:
        src = os.path.join(src_dir, file_name)
        print("copy '{}' -> '{}' directory".format(src, dst))
        shutil.copy2(src, dst)
