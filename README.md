HYDRO-BIDE
==========

Microbial ecological constraint-based community models

The repository contains:
* 1 Hydrodynamic BIDE (birth,immigration,death,emigration) simulation model based on neutral theory.
It contains aspects of conventional chemostate models (e.g. Monod, Droop).
* 1 Biogeographical simulation-based Monod model.
* 1 Biogeogrpahical simulation-based Droop (Cell Quota) Model.
* example files (.png) of figures the script creates
* files containing metagenomic data.  
  cow rumen data obtained from MG-RAST  http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=24 
* Module holding plotting and analysis functions

License
=======

HYDRO-BIDE is licensed under the open source MIT License

Copyright (c) 2013 Ken Locey

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
Contact information
-------------------
Ken Locey's email: ken@weecology.org and locey@biology.usu.edu

*Go to [Ken's website](http://kenlocey.weecology.org)*


Programs and packages used by the authors
-------------------------------

Python version 2.6.5 or higher
Numpy
Matplotlib
Cython (optional)

Cython
------
Most scripts will eventually have associated Cython setup files.


Cython is a language that makes writing C extensions for the Python language as easy as Python itself.
It is based on the well-known Pyrex, but supports more cutting edge functionality and optimizations.
*http://www.cython.org/*  
setup tutorial: *http://docs.cython.org/src/userguide/tutorial.html*

Most any Python script can be run at C speeds by:
* replacing the .py extension in your script with .pyx
* create a setup file. To save the user a hassle, we have made several setup files for scripts in this repository.
* run this command from within the folder containing the .pyx file and the setup file:


        python setup_thescriptsname.py build_ext --inplace

iPython and iPython Notebook
----------------------------

We are currently rescripting copies of repository files to be used in iPython Notebook. Running scripts in
iPython notebook will enable users to read comments in a convenient form, run blocks of code, and visualize
results between blocks of code. 
