import os
import shutil

# get list of notebooks
src_pth = os.path.join("..", "notebooks")
nb_files = [file_name
            for file_name in sorted(os.listdir(src_pth))
            if file_name.endswith(".ipynb")]

# create temporary directory
dst_pth = os.path.join("..", ".nb")
if os.path.isdir(dst_pth):
    shutil.rmtree(dst_pth)
os.makedirs(dst_pth)

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
if os.path.isdir(rtd_pth):
    print("deleting...'{}'".format(rtd_pth))
    shutil.rmtree(rtd_pth)

# copy the files
print("copying '{}' -> '{}'".format(dst_pth, rtd_pth))
shutil.copytree(dst_pth, rtd_pth)

# clean up temporary files
print("cleaning up...'{}'".format(dst_pth))
shutil.rmtree(dst_pth)