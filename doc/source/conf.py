# -*- coding: utf-8 -*-
# Py-ART documentation configuration file

import sys, os, re

# Check Sphinx version
import sphinx
if sphinx.__version__ < "1.0.1":
    raise RuntimeError("Sphinx 1.0.1 or newer required")

needs_sphinx = '1.0'

#----------------------------------------------------------------------------
# General configuration
#----------------------------------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.

sys.path.insert(0, os.path.abspath('../sphinxext'))

# Try to override the matplotlib configuration
try:
    import gen_rst
except:
    pass

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.mathjax',
              'sphinx.ext.autosummary', 'numpydoc']
# only include examples if the BUILD_PYART_EXAMPLES env. variable is set
if 'BUILD_PYART_EXAMPLES' in os.environ:
    extensions.append('gen_rst')

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'Py-ART'
copyright = u'2013-2019, Py-ART developers'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

import pyart
# The short X.Y version (including the .devXXXX suffix if present)
version = re.sub(r'^(\d+\.\d+)\.\d+(.*)', r'\1\2', pyart.__version__)
if 'dev' in version:
    # retain the .dev suffix, but clean it up
    version = re.sub(r'(\.dev\d*).*?$', r'\1', version)
else:
    # strip all other suffixes
    version = re.sub(r'^(\d+\.\d+).*?$', r'\1', version)
# The full version, including alpha/beta/rc tags.
release = pyart.__version__
# full Py-ART version in CI built docs
if 'CI' in os.environ and os.environ['CI'] == 'true':
    version = release

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_templates/*']

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------

# The style sheet to use for HTML and HTML Help pages. A file of that name
# must exist either in Sphinx' static/ path, or in one of the custom paths
# given in html_static_path.
html_theme = 'classic'
html_style = 'scipy.css'

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "Py-ART Documentation"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
html_file_suffix = '.html'

# Output file base name for HTML help builder.
htmlhelp_basename = 'pyart'

# ---------------------------------------------------------------------------
# LaTeX output
#----------------------------------------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'Py-ART.tex', u'Py-ART documentation',
   u'Py-ART developers', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True

#----------------------------------------------------------------------------
# Numpydoc extension
#----------------------------------------------------------------------------

# Numpy autodoc attributes
numpydoc_show_class_members = True

#----------------------------------------------------------------------------
# Autosummary
#----------------------------------------------------------------------------

if sphinx.__version__ >= "0.7":
    import glob
    #autosummary_generate = glob.glob("*.rst")
    autosummary_generate = True

#----------------------------------------------------------------------------
# Source code links
#----------------------------------------------------------------------------

# these functions borrowed from the scipy project
import inspect
from os.path import relpath, dirname
import pyart

for name in ['sphinx.ext.linkcode', 'linkcode', 'numpydoc.linkcode']:
    try:
        __import__(name)
        extensions.append(name)
        break
    except ImportError:
        pass
else:
    print "NOTE: linkcode extension not found -- no links to source generated"

def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
    """
    if domain != 'py':
        return None

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except:
            return None

    try:
        fn = inspect.getsourcefile(obj)
    except:
        fn = None
    if not fn:
        try:
            fn = inspect.getsourcefile(sys.modules[obj.__module__])
        except:
            fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.findsource(obj)
    except:
        lineno = None

    if lineno:
        linespec = "#L%d" % (lineno + 1)
    else:
        linespec = ""

    fn = relpath(fn, start=dirname(pyart.__file__))

    return "http://github.com/ARM-DOE/pyart/blob/master/pyart/%s%s" % (fn, linespec)
