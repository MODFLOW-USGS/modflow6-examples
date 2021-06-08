import os
import sys
import shutil
from subprocess import Popen, PIPE

zip = False
run = False
for arg in sys.argv:
    if arg.lower() in ("-z", "--zip"):
        zip = True
    elif arg.lower() in ("-r", "--run"):
        run = True

# clean examples directory if it exists
pth = os.path.join("..", "examples")
if os.path.isdir(pth):
    shutil.rmtree(pth)
os.makedirs(pth, exist_ok=True)

# run all of the scripts
pth = os.path.join("..", "scripts")
scripts = [
    file_name
    for file_name in sorted(os.listdir(pth))
    if file_name.endswith(".py") and file_name.startswith("ex-")
]
for s in scripts:
    args = ["python", s]
    if not run:
        args += ["--no_run", "--no_plot"]
    print(" ".join(args))
    proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print("Errors:\n{}".format(stderr.decode("utf-8")))

# zip the input files
if zip:
    print("zipping files in the 'examples' directory")
    sdir = os.getcwd()
    os.chdir("..")
    shutil.make_archive("modflow6-examples", "zip", "examples")
    os.chdir(sdir)
