import os
import sys
import shutil
import pytest

# path to notebooks
src_pth = os.path.join("..", "notebooks")
dst_pth = os.path.join("..", ".nb")
rtd_pth = os.path.join("..", ".doc", "_notebooks")

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
    nb_files = [
        file_name
        for file_name in sorted(os.listdir(src_pth))
        if file_name.endswith(".ipynb")
    ]


def clean_files():
    if os.path.isdir(dst_pth):
        shutil.rmtree(dst_pth)


@pytest.mark.parametrize(
    "file_name",
    nb_files,
)
def test_run_notebooks(file_name):
    # make paths if they do not exist
    for dir_path in (
        dst_pth,
        rtd_pth,
    ):
        os.makedirs(dir_path, exist_ok=True)

    # set src, dst, and rtd paths
    src = os.path.join(src_pth, file_name)
    dst = os.path.join(dst_pth, file_name)
    rtd = os.path.join(rtd_pth, file_name)

    # remove dst if it exists
    if os.path.isfile(dst):
        print(f"removing '{dst}'")
        os.remove(dst)

    arg = (
        "jupytext",
        "--to ipynb",
        "--from ipynb",
        "--execute",
        "-o",
        dst,
        src,
    )
    print(" ".join(arg))
    os.system(" ".join(arg))

    # remove rtd if it exists
    if os.path.isfile(rtd):
        print(f"removing '{rtd}'")
        os.remove(rtd)

    # copy dst to rtd
    print(f"copying '{dst}' -> '{rtd}'")
    shutil.copyfile(dst, rtd)


if __name__ == "__main__":
    clean_files()

    # clean up ReadtheDocs path
    if clean_sphinx:
        if os.path.isdir(rtd_pth):
            print(f"deleting...'{rtd_pth}'")
            shutil.rmtree(rtd_pth)

    for file_name in nb_files:
        test_run_notebooks(file_name)

    # clean up temporary files
    print(f"cleaning up...'{dst_pth}'")
    shutil.rmtree(dst_pth)
