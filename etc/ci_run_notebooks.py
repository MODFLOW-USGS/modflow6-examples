import os
import sys
import shutil

# path to notebooks
src_pth = os.path.join("..", "notebooks")

# parse command line arguments for notebook to create
nb_files = None
clean_sphinx = True
for idx, arg in enumerate(sys.argv):
    if arg in ("-f", "--file"):
        file_name = sys.argv[idx + 1]
        if not file_name.endswith(".ipynb"):
            file_name += ".ipynb"
        pth = os.path.join(src_pth, file_name)
        if os.path.isfile(pth):
            nb_files = [file_name]
            clean_sphinx = False

# get list of notebooks
if nb_files is None:
    nb_files = [file_name
                for file_name in sorted(os.listdir(src_pth))
                if file_name.endswith(".ipynb")]

# create temporary directory
dst_pth = os.path.join("..", ".nb")
if os.path.isdir(dst_pth):
    shutil.rmtree(dst_pth)
os.makedirs(dst_pth, exist_ok=True)

# run the notebooks
for file_name in nb_files:
    src = os.path.join(src_pth, file_name)
    dst = os.path.join(dst_pth, file_name)
    arg = ("jupytext",
           "--to ipynb",
           "--from ipynb",
           "--execute",
           "-o",
           dst,
           src)
    print(" ".join(arg))
    os.system(" ".join(arg))

# clean up ReadtheDocs path
rtd_pth = os.path.join("..", ".doc", "_notebooks")
if clean_sphinx:
    if os.path.isdir(rtd_pth):
        print("deleting...'{}'".format(rtd_pth))
        shutil.rmtree(rtd_pth)

# copy the files
if clean_sphinx:
    print("copying '{}' -> '{}'".format(dst_pth, rtd_pth))
    shutil.copytree(dst_pth, rtd_pth)
else:
    for file_name in nb_files:
        src = os.path.join(dst_pth, file_name)
        dst = os.path.join(rtd_pth, file_name)
        print("copying '{}' -> '{}'".format(src, dst))
        if os.path.isfile(dst):
            os.remove(dst)
        shutil.copyfile(src, dst)

# clean up temporary files
print("cleaning up...'{}'".format(dst_pth))
shutil.rmtree(dst_pth)