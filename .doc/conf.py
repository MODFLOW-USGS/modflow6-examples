# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import pymake

# -- processed notebook path -------------------------------------------------
nb_pth = os.path.join("_notebooks")

# -- download the completed notebooks -------------------------------------------------------
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    key = "modflow6-notebooks.zip"
    url = pymake.get_repo_assets("MODFLOW-USGS/modflow6-nightly-build")[key]
    pymake.download_and_unzip(url, nb_pth, verbose=True)

# -- get list of notebooks ---------------------------------------------------
nb_files = [os.path.join(nb_pth, file_name) for file_name in sorted(os.listdir(nb_pth))
            if file_name.endswith(".ipynb")]

# -- get list of python files ------------------------------------------------
pth = os.path.join("..", "scripts")
py_files = [os.path.join(pth, file_name) for file_name in sorted(os.listdir(pth)) if
            file_name.startswith("ex-") and file_name.endswith(".py")]

# -- process intro text out of python scripts --------------------------------
intro_text = []
for idx, fpth in enumerate(py_files):
    with open(fpth) as f:
        lines = f.read().splitlines()
    iend = 0
    for jdx, line in enumerate(lines):
        if not line.startswith("#"):
            iend = jdx
            break
    intro = []
    heading = 0
    for line in lines[:jdx + 1]:
        line = line[1:].strip()
        if len(line) < 1:
            continue
        elif line.startswith("##"):
            heading = 3
            line = line[2:]
        elif line.startswith("#"):
            heading = 2
            line = line[1:]
        else:
            if heading > 0:
                if heading == 1:
                    s = "="
                elif heading == 2:
                    s = "-"
                elif heading == 3:
                    s = "^"
                else:
                    s = " "
                intro.append(79 * "{}".format(s))
                intro.append("")
                intro.append("")
            heading = 0
        if len(line) > 0:
            intro.append(line)
    intro_text.append(intro)

# -- Build examples.rst for notebooks to .doc --------------------------------
f = open("examples.rst", "w")
for idx, fpth in enumerate(nb_files):
    file_name = os.path.basename(fpth)
    rst_pth = os.path.join("_notebooks", file_name)
    lines = ""
    for ex_list in intro_text[idx]:
        lines += "{}\n".format(ex_list.strip())
    lines += "Contents:\n\n"
    lines += "\n.. toctree::\n"
    lines += "   :maxdepth: 2\n\n"
    lines += "   {}\n\n\n".format(rst_pth)
    f.write(lines)
f.close()

# -- Project information -----------------------------------------------------

project = 'MODFLOW 6 Example Problems'
copyright = '2020, Langevin, C.D., Morway, E.D., and Hughes, J.D.'
author = 'Langevin, C.D., Morway, E.D., and Hughes, J.D.'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_console_highlighting",  # lowercase didn't work
    "sphinx.ext.autosectionlabel",
    "nbsphinx",
    "nbsphinx_link",
    "recommonmark",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
# The encoding of source files.
source_encoding = "utf-8"

# The master toctree document.
master_doc = "index"

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/MODFLOW-USGS/modflow6-examples",
    "use_edit_page_button": False,
}

autosummary_generate = True
numpydoc_show_class_members = False

html_context = {
    "github_user": "modflow6-examples",
    "github_repo": "modflow6-examples",
    "github_version": "master",
    "doc_path": ".doc",
}

# html_css_files = [
#     "css/custom.css",
# ]

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "modflow6-examples"
# html_favicon = ".._images/flopylogo.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".


# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = False

# Output file base name for HTML help builder.
htmlhelp_basename = "mf6exdoc"
